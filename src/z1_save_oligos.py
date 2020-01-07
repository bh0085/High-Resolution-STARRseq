#LOAD DATA FILES
from _config import OLIGO_LIBRARY, POSITIVE_CONTROLS_FILE
import numpy as np
import pandas as pd, re

MODERATE_EXPR_MU = 15
HIGH_EXPR_MU = 50
EXP_MODERATE_EXPR_MU = 3.3



oligos_lib = OLIGO_LIBRARY
oligos_lib["id"] = oligos_lib.index
oligos_lib["starts"] = oligos_lib.chromosome_info.apply(lambda x: int(re.compile("(\d+)-(\d+)").search(x).groups()[0]))
oligos_lib["starts"] = oligos_lib.starts-oligos_lib.starts.min()

candid_regex = re.compile("candid_(\d+)")
neg_regex = re.compile("Neg_(\d+)")
mutant_regex = re.compile("mut(\d+)$")

oligos = pd.read_csv("../out/0628_d1a_all_oligo_info.csv")
oligos_lib["locus_ids"] = oligos_lib.RefSeqID.apply(lambda x: int(candid_regex.search(x).groups()[0]) if candid_regex.search(x) else neg_regex.search(x).groups()[0] )
oligos_lib["mutant_num"] = oligos_lib.RefSeqID.apply(lambda x: int(mutant_regex.search(x).groups()[0]) if mutant_regex.search(x) else 0 )
oligos_lib["is_cand"] = oligos_lib.RefSeqID.apply(lambda x: True if candid_regex.search(x) else False )
oligos_lib["is_neg"] = oligos_lib.RefSeqID.apply(lambda x: True if neg_regex.search(x) else False )

unique_starts = pd.DataFrame(pd.Series(oligos_lib.starts.unique()).rename("starts"))
start_indexes = unique_starts.reset_index().set_index("starts")

oligos_lib["start_index"] = oligos_lib.starts.apply(lambda x: start_indexes.loc[x])



def save_oligos():
    #DEFINE CONSTANTS

    global oligos
    oligos = oligos.set_index("oligo")
    oligos = oligos.join(oligos_lib, on ="oligo")#.dropna()
    oligos["mu"] = oligos.n_transcripts / oligos.n_bcs
    oligos["start_index"] = oligos.starts.apply(lambda x: start_indexes.loc[x])
    chrom_info = oligos.groupby("chromosome_info").first().apply(lambda x: pd.Series(re.compile("(?P<chrom>[^:]*):(?P<gstart>\d+)-(?P<gend>\d+)").search(x.name).groupdict()),axis=1)
    oligos = oligos.reset_index().merge(chrom_info,on="chromosome_info").set_index("oligo")
    oligos["gstart"] = oligos.gstart.astype(np.int32)
    oligos["gend"] = oligos.gend.astype(np.int32)

    oligo_exp_info= pd.read_csv("../out/0628_d1a_all_exp_oligo_info.csv")
    oligos_by_exp = oligo_exp_info
    oligos_by_exp = oligos_by_exp.join(oligos_lib, on ="oligo")#.dropna()
    oligos_by_exp["mu"] = oligos_by_exp.n_transcripts / oligos_by_exp.n_bcs
    mean_mu = oligos_by_exp.mu.mean()
    mu_by_exp = oligos_by_exp.groupby("exp").mu.mean()

    oligos_by_exp["exp_mu"] =  oligos_by_exp.exp.apply(lambda x: mu_by_exp.loc[x])
    oligos_by_exp["exp_norm_mu"]= (oligos_by_exp["mu"] / oligos_by_exp["exp_mu"]) * mean_mu
    oligos_by_exp["exp_nm"] = oligos_by_exp.apply(lambda x:re.compile("_BR.").sub("",x.exp),axis=1)          
    oligos_by_exp["rep"] = oligos_by_exp.apply(lambda x:np.int32(re.compile("_BR(.)").search(x.exp).group(1)) if "BR" in x.exp else np.nan,axis=1)
    oligos_by_exp= oligos_by_exp.loc[oligos_by_exp.exp != "other"]

    oligos_by_exp["start_index"] = oligos_by_exp.starts.apply(lambda x: start_indexes.loc[x])
    chrom_info = oligos.groupby("chromosome_info").first().apply(lambda x: pd.Series(re.compile("(?P<chrom>[^:]*):(?P<gstart>\d+)-(?P<gend>\d+)").search(x.name).groupdict()),axis=1)

    oligos_by_exp = oligos_by_exp.reset_index().merge(chrom_info,on="chromosome_info").set_index(["exp","oligo"])
    oligos_by_exp["gstart"] = oligos_by_exp.gstart.astype(np.int32)
    oligos_by_exp["gend"] = oligos_by_exp.gend.astype(np.int32)

    oligos_by_exp = oligos_by_exp.reset_index()

    mut5s = oligos_by_exp.loc[lambda x: x.mutant_num == 0].copy()
    mut5s.mutant_num = 5
    mut5s.starts = mut5s.starts + 30
    oligo_map = mut5s.drop_duplicates("starts").set_index("starts").oligo + oligos_by_exp.oligo.max()
    mut5s["oligo"] = mut5s.starts.apply(lambda x: oligo_map.loc[x])
    oligos_by_exp= pd.concat([oligos_by_exp, mut5s])

    #processed data files
    oligos.to_csv("../out/0707_STARRSEQ_oligo_stats.csv")
    oligos_by_exp.to_csv("../out/0707_STARRSEQ_oligos_stats_by_experiment_all.csv")

def load_oligos():
    return [pd.read_csv("../out/0707_STARRSEQ_oligo_stats.csv",index_col="oligo"),
          pd.read_csv("../out/0707_STARRSEQ_oligos_stats_by_experiment_all.csv",index_col=["exp","oligo"])]

def load_oligos_plus():
    oligos,oligos_by_exp = load_oligos()    
    oligos_by_exp["norm_mu"] = oligos_by_exp["exp_norm_mu"]
    oligos_by_exp = oligos_by_exp.filter(regex='^(?!exp|\\.).*')
    oligos_by_exp["exp_nm"] = oligos_by_exp.index.get_level_values("exp").to_series().apply(lambda x:re.compile('(.*)_BR').search(x).groups()[0]).values
    oligos_by_exp = oligos_by_exp.reset_index()

    oligos_by_exp["exp_ct"] = oligos_by_exp.exp_nm.apply(lambda x:re.compile("(U2OS|DLD1|HCT116)").search(x).groups()[0])
    oligos_by_exp["exp_type"] = oligos_by_exp.exp_nm.apply(lambda x:"U2OS_NFKB" if "NFKB" in x 
                                                        else ("HCT116_GEM") if "Gem" in x
                                                        else re.compile("(U2OS|DLD1|HCT116)").search(x).groups()[0]+"_WT")

   
    oligos_by_exp["mutant_start"] = oligos_by_exp.reset_index().apply(lambda x: x.starts +x.mutant_num*30,axis=1)

    #oligos["mutant_start"] = oligos.mutant_start + 30
    #oligos_by_exp["mutant_start"] = oligos_by_exp.apply(lambda x: np.nan if x.mutant_num == 0 else x.mutant_start,axis=1)
    #oligos["mutant_start"] = oligos.apply(lambda x:np.nan if x.mutant_num == 0 else x.mutant_start,axis=1)

    ranks = oligos_by_exp.loc[lambda x:x.mutant_num == 0].groupby("exp_type").\
        apply(lambda x:x.groupby("oligo").mu.mean().sort_values(ascending=False).\
            reset_index().reset_index().set_index("oligo").\
            rename({"index":"ranksort"},axis="columns").\
            ranksort).unstack(level=0).rename(lambda x: x+"_rank",axis = "columns")
    oligos_by_exp = oligos_by_exp.join(ranks, on="oligo")

    return (None, oligos_by_exp)