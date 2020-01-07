import pandas as pd
import numpy as np
import re
from scipy.stats import ks_2samp
from scipy import stats

def get_ofs_opts(df):
    starts = df.index.get_level_values("starts")
    return np.array([i for i,start in enumerate(starts) if  (len(starts[i:i+4]) ==len(starts[i+1:i+5]) ) and (sum(starts[i:i+4] - starts[i+1:i+5]) == -120)])

def get_ofs_idxs(df):
    starts = df.index.get_level_values("starts")
    ofs_range = np.array([i for i,start in enumerate(starts) if  (len(starts[i:i+4]) ==len(starts[i+1:i+5]) ) and (sum(starts[i:i+4] - starts[i+1:i+5]) == -120)])
    return df.index.get_level_values("starts")[ofs_range+4]
    

def compute_filters(all_obe, hq = False):
    
    opts_df = all_obe.unstack(level=1).groupby(level=1).apply(
        lambda df: get_ofs_opts(df))
    idxs_df = all_obe.unstack(level=1).groupby(level=1).apply(
        lambda df: get_ofs_idxs(df))

    filters = pd.concat(
    [all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum((df.values[ofs:ofs+4,:] * (
    np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(1 - np.diag([1,1,1,1]))],
    axis=1))
        ))
    for ofs in get_ofs_opts(df)],index = idxs_df.loc[df.name])
    ).rename("allothers") / 16,
    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum((df.values[ofs:ofs+4,:] * (
    np.concatenate([
        (np.zeros((4,1))),
        np.fliplr( np.diag([1,1,1,1]))],
    axis=1))
        ))  
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("onlyablations") / 4,

    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum(
    (df.values[ofs:ofs+5,:][np.nonzero(
                1-np.array(
                    [[0,0,0,0,1],
                    [0,0,0,1,1],
                    [0,0,1,1,0],
                    [0,1,1,0,0],
                    [0,1,0,0,0]]))]
        ))
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("allothers2") / 17,
    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nanstd(df.values[ofs:ofs+4,:][np.nonzero (
    np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(1 - np.diag([1,1,1,1]))],
    axis=1))
                                        ])
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("allothersstd"),
    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nanstd(df.values[ofs:ofs+4,:])  
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("std"),
    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum(
    (df.values[ofs:ofs+5,:][np.nonzero(
                np.array([[0,0,0,0,1],
                [0,0,0,1,1],
                [0,0,1,1,0],
                [0,1,1,0,0],
                [0,1,0,0,0]]))]
        ))  
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("onlyablations2") / 8,
    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum((df.values[ofs:ofs+4,:] * (
    np.concatenate([
        (np.zeros((4,1))),
        np.fliplr(1 - np.diag([1,1,1,1]))],
    axis=1))
        ))  
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("othermutants") / 12,

    all_obe.unstack(level=1).groupby(level=1).apply(
    lambda df:
    pd.Series([np.nansum((df.values[ofs:ofs+4,:] * (
    np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(0* np.diag([1,1,1,1]))],
    axis=1))
        )) 
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])
    ).rename("onlywildtype") / 4,
    all_obe.unstack(level=1).groupby(level=1).apply(lambda df:

    pd.Series([ks_2samp( df.values[ofs:ofs+4,:][np.nonzero(   np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(1-np.diag([1,1,1,1]))],
    axis=1))] ,
            df.values[ofs:ofs+4,:][np.nonzero(  np.concatenate([
        (np.zeros((4,1))+0),
        np.fliplr(np.diag([1,1,1,1]))],
    axis=1))])[1]
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])).rename("ks_pval"),

    all_obe.unstack(level=1).groupby(level=1).apply(lambda df:

    pd.Series([ks_2samp( df.values[ofs:ofs+5,:][np.nonzero(
                np.array([[0,0,0,0,1],
                [0,0,0,1,1],
                [0,0,1,1,0],
                [0,1,1,0,0],
                [0,1,0,0,0]]))],
    df.values[ofs:ofs+5,:][np.nonzero(
                1- np.array([[0,0,0,0,1],
                [0,0,0,1,1],
                [0,0,1,1,0],
                [0,1,1,0,0],
                [0,1,0,0,0]]))]
    )[1]
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])).rename("ks2_pval"),


    all_obe.unstack(level=1).groupby(level=1).apply(lambda df:

    pd.Series([stats.ttest_ind( df.values[ofs:ofs+4,:][np.nonzero(   np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(1-np.diag([1,1,1,1]))],
    axis=1))] ,
            df.values[ofs:ofs+4,:][np.nonzero(  np.concatenate([
        (np.zeros((4,1))+0),
        np.fliplr(np.diag([1,1,1,1]))],
    axis=1))])
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])).apply(lambda x:pd.Series([x[0],x[1]])).rename({0:"ttest_ttval",1:"ttest_pval"},axis=1),

    all_obe.unstack(level=1).groupby(level=1).apply(lambda df:

    pd.Series([stats.ttest_ind( 
        df.values[ofs:ofs+5,:][np.nonzero(
                1- np.array([[0,0,0,0,1],
                [0,0,0,1,1],
                [0,0,1,1,0],
                [0,1,1,0,0],
                [0,1,0,0,0]]))],
        df.values[ofs:ofs+5,:][np.nonzero(
                np.array([[0,0,0,0,1],
                [0,0,0,1,1],
                [0,0,1,1,0],
                [0,1,1,0,0],
                [0,1,0,0,0]]))]
    )
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])).apply(lambda x:pd.Series([x[0],x[1]])).rename({0:"ttest2_ttval",1:"ttest2_pval"},axis=1),
    all_obe.unstack(level=1).groupby(level=1).apply(lambda df:

    pd.Series([ks_2samp( df.values[ofs:ofs+4,:][np.nonzero(   np.concatenate([
        (np.zeros((4,1))+1),
        np.fliplr(1-np.diag([1,1,1,1]))],
    axis=1))] ,
            df.values[ofs:ofs+4,:][np.nonzero(  np.concatenate([
        (np.zeros((4,1))+0),
        np.fliplr(np.diag([1,1,1,1]))],
    axis=1))])[0]
    for ofs in opts_df.loc[df.name]],index = idxs_df.loc[df.name])).rename("ks_stat")
    ],axis = 1)

    filters.index.names = ["analysis_group_key","starts"]
    filters["actual_starts"] = filters.index.get_level_values("starts") + 30
    filters = filters.reset_index(level=1, drop = True).set_index("actual_starts",append=True)
    filters.index.names = ["analysis_group_key","mutant_start_position"]

    filters["mutdiff"] = filters["othermutants"] - filters["onlyablations"]
    filters["wtdiff"] = filters["onlywildtype"] - filters["onlyablations"]\

    filters["othersdiff"] = filters["allothers"] - filters["onlyablations"]
    filters["rank_mutdiff"] = filters[["mutdiff"]].join(
            filters[["mutdiff"]].groupby("mutant_start_position").mean().reset_index().sort_values("mutdiff",ascending = False).reset_index().rename_axis("rank",axis="index").reset_index().set_index(["mutant_start_position"])["rank"],
        on="mutant_start_position")["rank"]
    filters["rank_ao"] = filters[["othersdiff"]].join(
            filters[["othersdiff"]].groupby("mutant_start_position").mean().reset_index().sort_values("othersdiff",ascending = False).reset_index().rename_axis("rank",axis="index").reset_index().set_index(["mutant_start_position"])["rank"],
        on="mutant_start_position")["rank"]
    filters["rank_ao_dld1"] = filters[["othersdiff"]].loc[lambda x: x.index.get_level_values(0).str.contains("DLD")].join(
            filters[["othersdiff"]].groupby("mutant_start_position").mean().reset_index().sort_values("othersdiff",ascending = False).reset_index().rename_axis("rank",axis="index").reset_index().set_index(["mutant_start_position"])["rank"],
        on="mutant_start_position")["rank"]
    filters["rank_ao_u2os"] = filters[["othersdiff"]].loc[lambda x: x.index.get_level_values(0).str.contains("U2OS")].join(
            filters[["othersdiff"]].groupby("mutant_start_position").mean().reset_index().sort_values("othersdiff",ascending = False).reset_index().rename_axis("rank",axis="index").reset_index().set_index(["mutant_start_position"])["rank"],
        on="mutant_start_position")["rank"]
    filters["rank_ao_hct116"] = filters[["othersdiff"]].loc[lambda x: x.index.get_level_values(0).str.contains("HCT")].join(
            filters[["othersdiff"]].groupby("mutant_start_position").mean().reset_index().sort_values("othersdiff",ascending = False).reset_index().rename_axis("rank",axis="index").reset_index().set_index(["mutant_start_position"])["rank"],
        on="mutant_start_position")["rank"]


    filters["filterchange"] = np.max(
        [(filters.onlyablations - filters.allothers),
        (filters.onlyablations2 - filters.allothers2)]
    )
    filters["ks_1or2"] = filters.apply(lambda x: (x.ks_pval < .05) | (x.ks2_pval < .05),axis =1)
    filters["log_ks_pval"] = np.log(filters.ks_pval)


    if hq:
        filters_hq = filters.loc[lambda x: x.ks_1or2]
        return filters_hq
    else:
        return filters


def compute_quantiles(filters):


    qvals = [.75,.8,.85,.9,.95,.98]
    quantiles = pd.DataFrame()
    quantile_counts = pd.DataFrame()

    for quantile in qvals:
        these_filters = filters.copy()
        
        wt_quantile = quantile
        ablation_quantile = quantile
        change_quantile = quantile

        #USES ALL EXPRESSION VALUES
        mubar_wt_cutoffs = filters.groupby("analysis_group_key").onlywildtype.quantile(wt_quantile)
        mubar_ablation_cutoffs = filters.groupby("analysis_group_key").onlyablations.quantile(ablation_quantile)

        these_filters["filterchange"] =  these_filters.apply(lambda x:
                                                        np.max(np.abs([x.onlyablations - x.allothers,
                                                                x.onlyablations2 - x.allothers2,])),axis=1)
        filter_change_cutoff = these_filters.filterchange.groupby("analysis_group_key").quantile(change_quantile)

        these_filters["ablation_mu_filtered"] = these_filters.apply(
            lambda y:(1 if 
                    (y.onlyablations > mubar_ablation_cutoffs.loc[y.name[0]]) &
                        (y.onlyablations > y.onlywildtype) else 0),axis=1)
        
        these_filters["wt_filtered"] = these_filters.apply(lambda y:(1 if
                    (y.onlywildtype > mubar_wt_cutoffs.loc[y.name[0]]) &
                        (y.onlywildtype > y.onlyablations) else 0),axis=1)
                                                                    
                                                                    
        these_filters["change_filtered"] = these_filters.apply(lambda y:(1 if (y.filterchange > filter_change_cutoff.loc[y.name[0]]) else 0),axis=1)
        these_filters["both_filtered"] = these_filters.wt_filtered * these_filters.change_filtered
        these_filters["both_ablation_filtered"] = these_filters.ablation_mu_filtered * these_filters.change_filtered

        these_filters["filter_color"] = these_filters.apply(
            lambda x:(6 if (x.both_filtered  and x.both_ablation_filtered)
                    else(5 if x.both_ablation_filtered
                            else(4 if x.ablation_mu_filtered
                                else(3 if x.both_filtered
                                    else(2 if x.change_filtered 
                                            else (1 if x.wt_filtered 
                                                else 0) ))))),axis=1)

        quantiles = quantiles.append(these_filters.assign(quantile = quantile))
        out = pd.Series([1,2,3]).apply(lambda n:(these_filters["filter_color"].loc[lambda x:x>0].unstack(level=0).fillna(0) == n).sum()).assign(quantile = quantile)
        quantile_counts =quantile_counts.append(out)
    return quantiles
