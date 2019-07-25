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
from _config import DATA_DIR, OUT_PLACE, N_SPLITS, QSUBS_DIR, OLIGO_LIBRARY, POSITIVE_CONTROLS_FILE,EXP_NAMES
import _config

#IO DIRECTORY CONFIG
NAME = util.get_fn(__file__)
OUT_DIR = os.path.join(OUT_PLACE, NAME)
util.ensure_dir_exists(OUT_DIR)

##
# CUSTOM CODE
##


#load oligo library from experimental design
oligo_lib = OLIGO_LIBRARY
oligo_lib["id"] = oligo_lib.index



def process_oligo(oligo_id_start,oligo_id_end):

    MERGED_INP_DIR = os.path.join(OUT_PLACE,"c1b_merge_tx_oligos_better")
    files = os.listdir(MERGED_INP_DIR)
    all_results2 = None

    for i,f in enumerate(files):
        #if f != "merged_56_7.csv": continue
        if "60" in f: 
            print( f"skipping {f}")
        #if i%25 == 0:
        print(f"{i} of {len(files)}")

        #print(i)
        #if i > 10:break
            
        INP_FILE = os.path.join(MERGED_INP_DIR,f)
        merged_results = pd.read_csv(INP_FILE,index_col=["oligo","bc","umi","exp"])
        print(merged_results.iloc[:10])
        merged_results.sort_index(level = 0, inplace=True)
        oligo_results = merged_results.loc[oligo_id_start:oligo_id_end]
        print("done loading")
        if all_results2 is None:
            all_results2 = oligo_results
        else:
            all_results2 = all_results2.add(oligo_results,fill_value =  0)
        print("done adding")
    all_results2.to_csv(os.path.join(OUT_DIR,f"merged_oligos_{oligo_id_start}.csv"))


##
# QSUB CODE
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
    print('Generating qsub scripts...')
    qsubs_dir = os.path.join(QSUBS_DIR ,NAME + '/')
    util.ensure_dir_exists(qsubs_dir)
    qsub_commands = []
    num_scripts = 0

    n_oligos_perbatch = 10

    for start in range(0,max(oligo_lib.index),n_oligos_perbatch):
        end = start + n_oligos_perbatch
        command = f'./{NAME}.py {start} {end}'
        script_id = NAME.split('_')[0]

        # Write shell scripts
        sh_fn = qsubs_dir + f'q_{script_id}_{start}.sh'
        with open(sh_fn, 'w') as f:
            f.write(f'#!/bin/bash\n{command}\n')
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
def main(oligo_id_start='', oligo_id_end=''):
  if oligo_id_start == '':
    gen_qsubs()
    return
  
  process_oligo(int(oligo_id_start),int(oligo_id_end))



if __name__ == '__main__':
  if len(sys.argv) > 1:
   main(*sys.argv[1:])
  else:
    main()    
