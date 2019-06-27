#!/cluster/bh0085/anaconda27/envs/py3/bin/python
##
# IMPORTS BOILERPLATE
##
import pandas as pd
import numpy as np
import itertools as it
import os,re, sys
from collections import defaultdict

sys.path.append("/cluster/bh0085")
from mybio import util
from _config import DATA_DIR, OUT_PLACE, N_SPLITS, QSUBS_DIR, EXP_NAMES, OLIGO_LIBRARY, POSITIVE_CONTROLS_FILE
import _config

#IO DIRECTORY CONFIG
NAME = util.get_fn(__file__)
OUT_DIR = os.path.join(OUT_PLACE, NAME)
util.ensure_dir_exists(OUT_DIR)

##
# CUSTOM CODE
##


#load oligo library from experimental design
oligos_lib = OLIGO_LIBRARY
oligos_lib["id"] = oligos_lib.index
positive_controls = pd.read_csv(POSITIVE_CONTROLS_FILE)


from _config import DATA_DIR, OUT_PLACE, N_SPLITS, QSUBS_DIR, EXP_NAMES
TX_INP_DIR = os.path.join(OUT_PLACE, "c0_bin_transcript_umis")
OLIGO_INP_DIR = os.path.join(OUT_PLACE, "b1_demultiplex_oligos")


#load oligo library from experimental design
oligos_lib = OLIGO_LIBRARY
oligos_lib["id"] = oligos_lib.index
positive_controls = pd.read_csv(POSITIVE_CONTROLS_FILE)




bc_substr_ids = dict([((s1+s2+s3) ,sum([     {"A":0,"T":1,"C":2,"G":3}[(s1+s2+s3)[i]]*(4**i) for i in range(3)])) for s1 in "ATGC" for s2 in "ATGC" for s3 in "ATGC"])
                     
#bc_id_substrs = dict([[(s1+s2+s3) ,sum([ {"A":0,"T":1,"C":2,"G":3}[(s1+s2+s3)[i]*(4**i)] for i in range(3)])] for s1 in "ATGC" for s2 in "ATGC" for s3 in "ATGC"])
bc_substr_inc = 4
bc_substr_max = len(bc_substr_ids)



def matchmaker(split, bc_substr_idx):
    
    bc_substr_idx = int(bc_substr_idx)

    merged_results = pd.DataFrame()
    tx_discovery_rates = []
    oligo_mapping_success_rates = []
    total_txs = 0
    total_oligos = 0
    #load barcode-matched umis from dictionary reads
    for oligo_file in os.listdir(OLIGO_INP_DIR):
        #file_format = "../out/b1_demultiplex_oligos/20190320-100XSTARRseq-Nextera-PE_HKNM2AFXY_S1_L001_3.json"
        if not ("_"+split+".json" in oligo_file):
            continue
        
        print(oligo_file)
        
        a_output2 = pd.read_csv(os.path.join(OLIGO_INP_DIR,oligo_file))
        a_output2 = a_output2.rename({"umi":"bc"}, axis = "columns")
        a_output_15 = a_output2.loc[a_output2.bc.str.len() == 15]
        total_oligos+= len(a_output2)
        
        cnt = 0

        for k, ogroup in a_output_15.groupby(a_output_15.bc.str.slice(0,3)):
            if "N" in k: 
              continue

            if (bc_substr_ids[k] < bc_substr_idx ) or ( bc_substr_ids[k] >= (bc_substr_idx + bc_substr_inc) ) :
              continue 
            cnt+=1

            for exp in EXP_NAMES:
        
                tx_file = os.path.join(TX_INP_DIR,"{}_{}.csv".format(exp,k))
                tx_subdict_df = pd.read_csv(tx_file)
                tx_subdict_df["exp"] = exp
                total_txs += len(tx_subdict_df)
                merged_results = merged_results.append(pd.merge(tx_subdict_df,ogroup,on="bc"))
                                        

                oligo_mapping_success_rates.append(float(len(merged_results)) / len(ogroup))
                tx_discovery_rates.append(float(len(merged_results)) / len(tx_subdict_df))


    merged_results.to_csv(os.path.join(OUT_DIR,f"merged_{bc_substr_idx}_{split}.csv"),index=False)




##
# QSUB CODE
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
    print('Generating qsub scripts...')
    qsubs_dir = QSUBS_DIR + NAME + '/'
    util.ensure_dir_exists(qsubs_dir)
    qsub_commands = []
    num_scripts = 0

    for split in range(N_SPLITS):
        for bc_substr_idx in range(0,bc_substr_max, bc_substr_inc):
            command = './%s.py %s %s' % (NAME, split, bc_substr_idx)
            script_id = NAME.split('_')[0]

            # Write shell scripts
            sh_fn = qsubs_dir + f'q_{script_id}_{bc_substr_idx}_{split}.sh'
            with open(sh_fn, 'w') as f:
                f.write('#!/bin/bash\n%s\n' % (command))
                num_scripts += 1
        
            # Write qsub commands
            qsub_commands.append(f'qsub -m e -wd {_config.SRC_DIR} {sh_fn}')
        
   # Save commands
    with open(qsubs_dir + '_commands.txt', 'w') as f:
        f.write('\n'.join(qsub_commands))

    print(f'Wrote {num_scripts} shell scripts to {qsubs_dir}')
    return


##
# Main
##
@util.time_dec
def main(split='', bc_substr_idx='' ):
  if split == '':
    gen_qsubs()
    return

  matchmaker(split, bc_substr_idx)
  return OUT_DIR


if __name__ == '__main__':
  if len(sys.argv) > 1:
   main(*sys.argv[1:])
  else:
    main()    
