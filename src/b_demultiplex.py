# 
from __future__ import division
import _config
from _config import OUT_DIR, EXP_DESIGN_2955
import sys, os, fnmatch, datetime, subprocess, imp
import numpy as np
from collections import defaultdict

import pandas as pd

import itertools as it
import re

# Default params
inp_dir = OUT_DIR + 'a2_split_3202/'
NAME = util.get_fn(__file__)
out_dir = os.path.join(_config.OUT_DIR, NAME)
util.ensure_dir_exists(out_dir)

exp_design=EXP_DESIGN_2955


exp_test_strs = {}
for k,r in exp_design.iterrows():
    exp_test_strs[r.Name] = "{0}+{1}".format(r["i7 Index reads"].strip(),r["i5 Index reads"].strip())

##
# Functions
##
def match(r1,r2, h1,h2):   
  for k, v in exp_test_strs.items():
    try:
        idx = h1.index(v)
        return k,r1
    except ValueError, e:
        continue
  return "other",r1

##
# primary
##
def demultiplex(split):
  #for inp_fn in [inp_dir + 'Undetermined_AH3W5GBGX9_S0_L00{0}_R1_001_{1}.fastq'.format(k2, split) for k in range(1,5)]:
  for name in list(exp_design["Name"]) + ['other']:
        util.ensure_dir_exists(out_dir + name)
        util.exists_empty_fn(out_dir + name + '/R1_%s.fa' % (split))
        util.exists_empty_fn(out_dir + name + '/R2_%s.fa' % (split))

  for snum, sgroup in it.groupby(
          sorted(os.listdir(inp_dir),key=lambda x:re.compile("(\d+)\.fastq").search(x).groups()[0])
      , key=lambda x:re.compile("(\d+)\.fastq").search(x).groups()[0]):
    
    if snum != split: continue
    for lnum,lgroup in it.groupby(
        sorted(list(sgroup), key=lambda x:int(re.compile("_L(\d+)").search(x).group(1)))
        , key=lambda x:int(re.compile("_L(\d+)").search(x).group(1))):
    

        
        fns = list(lgroup)
        print "LANE: {0}, FILES: {1}".format(lnum, fns)
        #print fns
        #print lnum
        read_files =dict([[ int(re.compile("R(\d+)").search(e).group(1)),e] for e in fns])

        inp_fn1 = os.path.join(inp_dir, read_files[1])
        inp_fn2 = os.path.join(inp_dir, read_files[2])
    
        lc = util.line_count(inp_fn1)
        num_bad_q, num_tot, num_other = 0, 0, 0
        timer = util.Timer(total = lc)
        i = -1

        
        with open(inp_fn1) as f1:
          with open(inp_fn2) as f2:
            while 1:
              i+=1
              if i % 1000000 ==0 : print "{0} records, ({1}%) [{2} bad] [{3} other]".format(i/4 , 100*float(i) / lc, num_bad_q,num_other)
              try: 
                line1 = f1.next()
                line2 = f2.next()
              except StopIteration, e:
                break
            
              if i % 4 == 0:
                h1 = line1.strip()
                h2 = line2.strip()
              if i % 4 == 1:
                r1 = line1.strip()
                r2 = line2.strip()
              if i % 4 == 3:
                num_tot += 1
                qs1 = line1.strip()
                qs2 = line2.strip()
                for qs in [qs1,qs2]:
                    quals = [ord(s)-33 for s in qs]
                    if np.mean(quals) < 30:
                      num_bad_q += 1
                      continue
    
                demultiplex_id, trimmed_read = match(r1,r2, h1,h2)
    
                #raise Exception()
            
                #print demultiplex_id
                if demultiplex_id == 'other':
                  num_other += 1
              
                #break
             
                out1_fn = out_dir +  '%s/R1_%s.fa' % (demultiplex_id, split)
                with open(out1_fn, 'a') as f:
                  f.write('>' + h1[1:] + '\n' + r1 + '\n')
                  
                out2_fn = out_dir +  '%s/R2_%s.fa' % (demultiplex_id, split)
                with open(out2_fn, 'a') as f:
                  f.write('>' + h2[1:] + '\n' + r2 + '\n')
          
              #timer.update()
    
        print 'Rejected %s fraction of reads' % (num_bad_q / num_tot)
  return

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print 'Generating qsub scripts...'
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  for idx in range(0, 15):
    command = 'python %s.py %s' % (NAME, idx)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, idx)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -m e -wd %s %s' % (_config.SRC_DIR, sh_fn))

  # Save commands
  with open(qsubs_dir + '_commands.txt', 'w') as f:
    f.write('\n'.join(qsub_commands))

  print 'Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir)
  return



##
# Main
##
@util.time_dec
def main(split = ''):
  if split == '':
    gen_qsubs()
    return

  demultiplex(split)
  return out_dir




if __name__ == '__main__':
  if len(sys.argv) > 1:
   main(sys.argv[1])
  else:
    main()
