import pandas as pd
import numpy as np


#LOADS BIOLOGICAL MOTIFS AND SCANS ALL SUBREGIONS FOR OCCURENCES
from pyfaidx import Fasta
sequences_fa = Fasta('/Users/ben/genomes/GRCh38.primary_assembly.genome.fa')
chrseq = str(sequences_fa["chr22"])
region_bounds=[ 38699734, 39291007]
background = dict([[l,chrseq[region_bounds[0]:region_bounds[1]].count(l) / len(chrseq[region_bounds[0]:region_bounds[1]])] for l in "ATGC"])



#load jaspar motifs
def save_jaspar():
    #load motifs
    from Bio.motifs import jaspar as mjaspar
    import Bio.motifs
    from Bio import motifs as bmotifs
    with open("../data/jaspar.pfm") as handle:
        jaspar_motifs = mjaspar._read_jaspar(handle)

    for j in jaspar_motifs:
         j.pseudocounts =1   

    print("done loading motifs")

    #identify motif clusters
    from sklearn.cluster import SpectralClustering, MeanShift, KMeans
    import numpy as np

    print("computing thresholds")
    
    jaspar = pd.DataFrame([pd.Series({"consensus":str(j.consensus),"name" :j.name,"pwm":j.pwm,"pssm":j.pssm}) for i,j in enumerate(jaspar_motifs)])
    print ("computing all thresholds") 
    #jaspar["threshold_balanced"] = jaspar.pssm.iloc[:10].apply(lambda x:x.distribution(background=background, precision=10**3).threshold_balanced(1))
    #print("computing fdr threshold")
    jaspar["threshold_fdr_005"] = None
    jaspar["threshold_fdr_05"] = None
    jaspar["threshold_balanced"] = None
    jaspar["threshold_patser"] = None
    for k,r in jaspar.iterrows():
        print(k)
        d =  r.pssm.distribution(background=background, precision=10**3)
        jaspar.at[k,"threshold_fdr_005"] =d.threshold_fpr(.005)
        jaspar.at[k,"threshold_fdr_05"] =d.threshold_fpr(.05)
        jaspar.at[k,"threshold_bal"] =d.threshold_balanced()
        jaspar.at[k,"threshold_patser"] =d.threshold_patser()

    #jaspar["threshold_fdr_005"] = jaspar.pssm.apply(lambda x:x.distribution(background=background, precision=10**3).threshold_fpr(.005))
    
    print("initializing clustering")
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
        print("clustering length ", l)
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

    print("cleaning up")
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

    print(f"""background model for region {region_bounds},:
    {background}""")

    jaspar["pssm_score_mean"] = jaspar.pssm.apply(lambda x: x.mean(background))
    jaspar["pssm_score_std"] = jaspar.pssm.apply(lambda x: x.std(background))
    jaspar["pssm_score_threshold"] = jaspar.pssm_score_mean + 2 * jaspar.pssm_score_std



    print("saving")
    jaspar.index.name="jaspar_id"
    jaspar.to_csv("../out/0708_jaspar.csv")



def load_jaspar():

    #load motifs
    from Bio.motifs import jaspar as mjaspar
    import Bio.motifs
    from Bio import motifs as bmotifs
    with open("../data/jaspar.pfm") as handle:
        jaspar_motifs = mjaspar._read_jaspar(handle)

    for j in jaspar_motifs:
         j.pseudocounts =1   

    jaspar = pd.read_csv("../out/0708_jaspar.csv", index_col = "jaspar_id")
    jaspar["pssm"] = pd.Series([j.pssm for j in jaspar_motifs])
    return jaspar