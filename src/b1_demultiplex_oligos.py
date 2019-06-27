import os, re, sys, json
from collections import defaultdict
import numpy as np, pandas as pd
import _config

from _config import  OLIGO_LIBRARY, A_BC_OLIGOS_OUT, OUT_PLACE as OUT_DIR
sys.path.append('/cluster/bh0085/')

from mybio import util

# Default params
inp_dir = OUT_DIR + 'a2_split_all/'
NAME = util.get_fn(__file__)
out_dir = os.path.join(_config.OUT_PLACE, NAME)
util.ensure_dir_exists(out_dir)


##
# Locality sensitive hashing
##


global names_targets 
names_targets = dict([(k,v.Sequences) for k,v in OLIGO_LIBRARY.iterrows()])

def build_targets_better_lsh():
  global names_targets
  lsh_dict = defaultdict(list)
  for nm in names_targets:
    target = names_targets[nm]
    kmers = get_lsh_kmers(target)
    for kmer in kmers:
      lsh_dict[kmer].append(nm)
  return lsh_dict

def get_lsh_kmers(target):
  kmer_len = 10
  kmers = []
  for idx in range(len(target) - kmer_len):
    kmer = target[idx : idx + kmer_len]
    kmers.append(kmer)
  return kmers

def find_best_designed_target(read, lsh_dict):
  kmers = get_lsh_kmers(read)
  scores = dict()
  for kmer in kmers:
    for exp in lsh_dict[kmer]:
      if exp not in scores:
        scores[exp] = 0
      scores[exp] += 1

  if len(scores) == 0:
    return []

  sorted_scores = sorted(scores, key = scores.get, reverse = True)
  best_score = scores[sorted_scores[0]]
  secondbest_score = scores[sorted_scores[1]]
  cand_idxs = sorted_scores[0]
  return cand_idxs, best_score,secondbest_score


#print difflib.get_close_matches(oligo_seq, OLIGO_LIBRARY.Sequences.values[:100] , n=1, cutoff=.01)

def featurize0(s):
    return np.array([s.count(l) for l in "ATGC"],dtype=np.float)

def featurize(s):
    return np.array([{"A":0.0,"T":1.0,"G":2.0,"C":3.0,"N":4.0}[l] for l in s])

def featurize_dict(s):
    return dict([(l,s.count(l)) for l in "ATGC"])
def selectbest(query, subjects, subject_features):
    qcounts0 = featurize(query)
    norm = np.sqrt(np.sum(qcounts0**2))
    qcounts = qcounts0 / norm
    dotprod = np.sum(qcounts*subject_features,1)
    second ,cand = np.argsort(dotprod)[-2:]
    
    return cand, dotprod[cand], dotprod[second]



def demultiplex( split,nm):

   lsh_dict = build_targets_better_lsh()
   
   
   num_tot = 0
   num_bad_q = 0
   matched_reads = []
   qualities = []
   notfound = 0
   num_valerr = 0
   
   bcs_by_oligo = dict()
   poorly_aligned = 0   

   subject_features0 = np.array([featurize(s) for s in  OLIGO_LIBRARY.Sequences.values]) / len(OLIGO_LIBRARY.Sequences[0])
   subject_features = subject_features0 / np.sqrt(np.sum(subject_features0**2,1))[:,np.newaxis]
   subjects = OLIGO_LIBRARY.Sequences.values
   
   bc_oligo_pairs = set({})
   flush_number = 0
   
   bc_oligo_counts = {}
   i = -1

   fn1 = nm+"_R1_001_{}.fastq".format(split)
   fn2 = nm+"_R2_001_{}.fastq".format(split)

   print(fn1)
   print(fn2)
   
   r1fn =os.path.join(inp_dir,fn1)
   r2fn =os.path.join(inp_dir,fn2)
   with open(r1fn) as f1:
       with open(r2fn) as f2:
           while 1:
               if i %10000 == 0: print(float(i) / 10000)
               
               i+=1
               try:
                   #for j in range(1+4*1000):
                   l1 = next(f1)
                   l2 = next(f2)
               except StopIteration:
                   break
   
   
               if i % 4 == 0:
                   h1 = l1.strip()
                   h2 = l2.strip()
               if i % 4 == 1:
                   r1 = l1.strip()
                   r2 = l2.strip()
               if i % 4 == 3:
                   num_tot += 1
                   qs1 =l1.strip()
                   qs2 =l2.strip()
                   
                   #...RIGHT NOW, ONLY CHECK R2 FOR QUALITY... R1 TENDS TO STINK. TOO LONG
                   quals = [ord(s)-33 for s in qs2]
                   if np.mean(quals) < 30:
                       num_bad_q += 1
                       continue
                     
                   try:
                       oligo_offset = r1.index("TGCACCGG")
                   except ValueError as v:
                       notfound+=1
                       continue
                   
                   oligo_start = oligo_offset + len("TGCACCGG")
                   oligo_seq = r1[oligo_start:oligo_start+150]
                   
                   r2_consensus = "AATTCGTCGA"
                   try:
                       bc_offset = r2.index(r2_consensus)
                   except ValueError:
                       notfound+=1
                       continue
                   bc_start = bc_offset + len(r2_consensus)
                   bc_seq = r2[bc_start:bc_start + 15]

                   try:
                       best,score2,secondbest_score2 = find_best_designed_target(oligo_seq,lsh_dict)
                   except ValueError:
                       num_valerr+=1
                       continue
                   
                   if score2 < 100 or ( score2-secondbest_score2) <10:
                       poorly_aligned += 1
                   else:
                       bc_oligo_counts[(bc_seq,best)]= bc_oligo_counts.get((bc_seq,best),0)+1
                       
   df = pd.DataFrame({"bc":k[0],"oligo":k[1],"count":v} for k,v in bc_oligo_counts.items())
   df.to_csv(os.path.join(out_dir,"{}_{}.json".format(nm,split)),index=False)


##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
    print('Generating qsub scripts...')
    qsubs_dir = _config.QSUBS_DIR + NAME + '/'
    util.ensure_dir_exists(qsubs_dir)
    qsub_commands = []

    num_scripts = 0
    for fn in os.listdir(inp_dir):
       basename,rnum,snum = re.compile("(.*)_R(\d)_001_(\d+).fastq").search(fn).groups()
       if int(rnum)==1:
         command = '/cluster/bh0085/anaconda27/envs/py3/bin/python %s.py %s %s' % (NAME, snum, basename)
         script_id = NAME.split('_')[0]

         # Write shell scripts
         sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id,basename,snum)
         with open(sh_fn, 'w') as f:
           f.write('#!/bin/bash\n%s\n' % (command))
         num_scripts += 1
     
         # Write qsub commands
         qsub_commands.append('qsub -m e -wd %s %s' % (_config.SRC_DIR, sh_fn))
     
   # Save commands
    with open(qsubs_dir + '_commands.txt', 'w') as f:
     f.write('\n'.join(qsub_commands))

    print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
    return



##
# Main
##
@util.time_dec
def main(split = '',nm =''):
  if split == '':
    gen_qsubs()
    return

  demultiplex(split,nm)
  return out_dir




if __name__ == '__main__':
  if len(sys.argv) > 1:
   main(sys.argv[1], sys.argv[2])
  else:
    main()
