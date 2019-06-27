# 
import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np
sys.path.append('/cluster/bh0085/')
from mybio import util
import pandas as pd


# Default params
inp_dirs = [_config.SHE3202_DIR, _config.SHE3439_DIR]
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
print(out_dir)
util.ensure_dir_exists(out_dir)

##
# Functions
##


def split(inp_fn, out_nm):
  inp_fn_numlines = util.line_count(inp_fn)

  num_splits = 30
  split_size = int(inp_fn_numlines / num_splits)
  if num_splits * split_size < inp_fn_numlines:
    split_size += 1
  while split_size % 4 != 0:
    split_size += 1

  
  print('Using split size %s' % (split_size))

  split_num = 0
  for idx in range(1, inp_fn_numlines, split_size):
    start = idx
    end = start + split_size  
    out_fn = out_dir + out_nm + '_%s.fastq' % (split_num)
    command = 'tail -n +%s %s | head -n %s > %s' % (  start,inp_fn, end - start, out_fn)
    split_num += 1
    print(command)

  return


##
# Main
##
@util.time_dec
def main():
  print(NAME)  
  
  # Function calls
  for inp_dir in inp_dirs:
    for fn in os.listdir(inp_dir):
      if fn[-5:] == 'fastq':
        split(inp_dir + fn, fn.replace('.fastq', ''))

  return


if __name__ == '__main__':
  main()
