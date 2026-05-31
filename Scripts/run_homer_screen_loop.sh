#!/bin/bash

# Run homer findMotifs.pl
findMotifs.pl ../Data/loop_sequence/wnt5a_loop.fa fasta homer_loop/ -find $HOME/homer/data/knownTFs/vertebrates/known.motifs > found_motifs_loop.txt

# Select only high score hits
awk 'NR==1 || ($NF + 0) > 8' found_motifs_loop.txt > high_confidence_motif_hits_loop.txt

# $NF - last field
# +0 - converts the field to numeric
