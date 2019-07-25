import pandas as pd
import numpy as np

#load jaspar motifs
def load_motifs():
    #load motifs
    from Bio.motifs import jaspar as mjaspar
    import Bio.motifs
    from Bio import motifs as bmotifs
    with open("../data/jaspar.pfm") as handle:
        jaspar_motifs = mjaspar._read_jaspar(handle)


    #identify motif clusters
    from sklearn.cluster import SpectralClustering, MeanShift, KMeans
    import numpy as np
    jaspar = pd.DataFrame([pd.Series({"consensus":str(j.consensus),"name" :j.name,"pwm":j.pwm,"pssm":j.pssm}) for i,j in enumerate(jaspar_motifs)])
    jaspar["spec2_cluster"] = -1
    jaspar["spec3_cluster"] = -1
    jaspar["ms_cluster"] = -1
    jaspar["km_cluster"] = -1
    jaspar["ms_center"] = None
    jaspar["km_center"] = None
    jaspar["km_cluster_name"] = None
    jaspar["spec2_cluster_name"] = None
    jaspar["spec3_cluster_name"] = None
    jaspar["ms_cluster_name"] = None

    #is this actually helpful?
    for l,g in jaspar.groupby(jaspar.consensus.str.len()):
        X = np.array(list(g.pwm.apply(lambda x: np.array([e for v in x.values() for e in v]))))

        #use meanshift clustering to find an ideal n_clusters
        ms_clustering = MeanShift().fit(X)
        jaspar.loc[g.index,"ms_cluster"] =[int(e) for e in ms_clustering.labels_]
        ms_centers = ms_clustering.cluster_centers_
        jaspar.loc[g.index,"ms_center"] = jaspar.loc[g.index]["ms_cluster"].apply(lambda x:ms_centers[x])

        #run spectral clustering (not used)
        spec2 = SpectralClustering(
            n_clusters = len(ms_centers),
            assign_labels="discretize",
            random_state=0).fit(X)
        jaspar.loc[g.index,"spec2_cluster"]  =[int(e) for e in spec2.labels_]

        
        #run spectral clustering (not used)
        spec3 = SpectralClustering(
            n_clusters = 8,
            assign_labels="discretize",
            random_state=0).fit(X)
        jaspar.loc[g.index,"spec3_cluster"]  =[int(e) for e in spec3.labels_]
        
        
        #initializing from meanshift, run kmeans clustering to find clusters and canonical motifs
        km_cluster = KMeans(
            n_clusters = len(ms_centers),
            random_state=0).fit(X)
        jaspar.loc[g.index,"km_cluster"] =[int(e) for e in km_cluster.labels_]
        km_centers = km_cluster.cluster_centers_
        jaspar.loc[g.index,"km_center"] = jaspar.loc[g.index]["km_cluster"].apply(lambda x:km_centers[x])

        closest = [ np.argmin([ np.sum(np.square(X[i] -cent))   for i in range(len(X))  ]) for cent in km_centers]
        jaspar.loc[g.index,"km_centroid_jaspar_id"] = jaspar.loc[g.index].km_cluster.apply(lambda i:g.iloc[ closest[int(i)] ].name)

    jaspar["len"] = jaspar.consensus.str.len()
    jaspar["spec2_cluster_name"] = jaspar.apply(lambda x:f"{x.len}_{x.spec2_cluster}",axis=1)
    jaspar["km_cluster_name"] = jaspar.apply(lambda x:f"{x.len}_{x.km_cluster}",axis=1)
    jaspar["spec3_cluster_name"] = jaspar.apply(lambda x:f"{x.len}_{x.spec3_cluster}",axis=1)
    jaspar["ms_cluster_name"] = jaspar.apply(lambda x:f"{x.len}_{x.ms_cluster}",axis=1)


    jaspar_km_names_to_ids =\
        jaspar.km_cluster_name.drop_duplicates().reset_index(drop=True)\
        .rename("cluster_string_name").reset_index().set_index("cluster_string_name")["index"].rename("cluster_id")
    jaspar_spec2_names_to_ids =\
        jaspar.spec2_cluster_name.drop_duplicates().reset_index(drop=True)\
        .rename("cluster_string_name").reset_index().set_index("cluster_string_name")["index"].rename("cluster_id")
    jaspar_spec3_names_to_ids =\
        jaspar.spec3_cluster_name.drop_duplicates().reset_index(drop=True)\
        .rename("cluster_string_name").reset_index().set_index("cluster_string_name")["index"].rename("cluster_id")
    jaspar_ms_names_to_ids =\
        jaspar.ms_cluster_name.drop_duplicates().reset_index(drop=True)\
        .rename("cluster_string_name").reset_index().set_index("cluster_string_name")["index"].rename("cluster_id")

    jaspar["km_cluster_id"] = jaspar. km_cluster_name.apply(lambda x: jaspar_km_names_to_ids.get(x))
    jaspar["spec2_cluster_id"] = jaspar. spec2_cluster_name.apply(lambda x: jaspar_spec2_names_to_ids.get(x))
    jaspar["spec3_cluster_id"] = jaspar.spec3_cluster_name.apply(lambda x: jaspar_spec3_names_to_ids.get(x))
    jaspar["ms_cluster_id"] = jaspar.ms_cluster_name.apply(lambda x: jaspar_ms_names_to_ids.get(x))
    jaspar.to_csv("../out/0708_jaspar.csv")


def run():
    jaspar = load_motifs()
    jaspar.to_csv("../out/0708_jaspar.csv")

def get_output_df():
    return pd.read_csv("../out/0708_jaspar.csv")
