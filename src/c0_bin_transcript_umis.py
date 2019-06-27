import pandas as pd
import numpy as np
import itertools as it
import os,re, sys
from collections import defaultdict

sys.path.append("/cluster/bh0085")
from mybio import util
from _config import DATA_DIR, OUT_PLACE, N_SPLITS, QSUBS_DIR, EXP_NAMES
import _config

#IO DIRECTORY CONFIG
NAME = util.get_fn(__file__)
out_dir = os.path.join(OUT_PLACE, NAME)
util.ensure_dir_exists(out_dir)
inp_dir = os.path.join(OUT_PLACE, "b0_demultiplex_transcripts")


#THIS FUNCTION GETS TRANSCRIPTS!!!!
def wfun(r,ds,fnames):
    out = []
    for f in fnames:
        if f[-2:]=="fa":
            out.append(os.path.join(r,f))
    return out
fa_found = [e for r,d,fs in os.walk("../out/b0_demultiplex_transcripts/",wfun) for e in wfun(r,d,fs) ]



def list_umis(pool_prefix, split,exp):
    
    f_infos = dict([[fn,re.compile("(?P<pool_prefix>[^/]*)/(?P<exp>[^/]*)/R(?P<rnum>[12])_(?P<split>\d+).fa").search(fn).groupdict()] 
                     for fn in fa_found])

    f1 = [k for k,v in f_infos.items() if v["rnum"] == "1" and v["split"]==split and v["exp"]==exp and v["pool_prefix"]==pool_prefix][0] 
    f2 = [k for k,v in f_infos.items() if v["rnum"] == "2" and v["split"]==split and v["exp"]==exp and v["pool_prefix"]==pool_prefix][0] 

    skipcount = 0
    nl = 0
    nacc = 0
    nnreject = 0
    umis = []
    t_umis = []

    print(f1)
    with open(f1) as f1open:
      with open(f2) as f2open:
        i = -1
        while 1:
            i+=1
            try:
                l1 = next(f1open)
                l2 = next(f2open)
            except StopIteration:
                break

            if i %2 == 0: continue
            if "N" in l2[:16]:
                #print(l2)
                #print("rejecting")
                nnreject +=1
                continue
            #print("ccepting")

            if i % 10001 == 0:
                print(i)
            umi = l2[:15].strip()
            umis.append(umi)          
            t_umi = l1[:10].strip()
            t_umis.append(t_umi)

            if( ">" in t_umi) or (">" in umi):
                raise Exception()
            nacc+=1
    return umis, t_umis, {"nacc":nacc, "nnreject":nnreject}


def run_splits( exp):
    pool_prefixes = os.listdir(inp_dir)
    all_bcs = pd.DataFrame()
    for pool_prefix in pool_prefixes:
        splits = range(N_SPLITS)
        for split in splits:
            total_count = 0 
            bcs,umis,stats = list_umis(pool_prefix,str(split) , exp)
            if len(bcs) > 0:
                all_bcs = all_bcs.append([{"bc":bc,"umi":umis[i],"exp":exp} for i,bc in enumerate(bcs)],ignore_index=True)
            else:
                print (f"EMPTY BCS FILE {pool_prefix}")

    print("GROUPING")
    for k,grp in all_bcs.groupby(all_bcs.bc.str.slice(0,3)):
        out_file = os.path.join(out_dir, f"{exp}_{k}.csv")
        grp.to_csv(out_file)

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
    print('Generating qsub scripts...')
    qsubs_dir = QSUBS_DIR + NAME + '/'
    util.ensure_dir_exists(qsubs_dir)
    qsub_commands = []

    num_scripts = 0

    for exp in EXP_NAMES:
            command = '/cluster/bh0085/anaconda27/envs/py3/bin/python %s.py %s' % (NAME, exp)
            script_id = NAME.split('_')[0]

            # Write shell scripts
            sh_fn = qsubs_dir + f'q_{script_id}_{exp}.sh'
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
def main(exp='' ):
  if exp == '':
    gen_qsubs()
    return

  run_splits( exp)
  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
   main(*sys.argv[1:])
  else:
    main()    
