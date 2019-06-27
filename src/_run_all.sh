
source activate py3

#SPLIT TRANSCRIPTS LIBRARY
python a_split.py  > a_split_cmds.txt
source a_split_cmds.txt

#SPLIT OLIGOS LIBRARY
python a2_split_all.py > a2_split_all_cmds.txt
source a2_split_all_cmds.txt

#DEMULTIPLEX TRANSCRIPTS
python b0_demultiplex_transcripts.py
source ../qsubs/b0_demultiplex_transcripts/_commands.txt

#DEMULTIPLEX TRANSCRIPTS
python b1_demultiplex_oligos.py
source ../qsubs/b1_demultiplex_transcripts/_commands.txt

#BIN TRANSCRIPTS INTO FILES BY BARCODE PREFIX
python c0_bin_transcript_umis.py
source ../qsubst/c0_bin_transcript_umis/_commands.txt

#MERGE OLIGOS AND TRANSCRIPTS BY PREFIX
./c1_merge_tx_oligos.py
source ../qsubs/c1_merge_tx_oligos/_commands.txt