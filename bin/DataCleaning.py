#!/usr/bin/env python3
#Author: Pulkit Anupam Srivastava
#Co-authors: Eric L. Hegg, Michigan State University
#	     Brian G. Fox, University of Wisconsin-Madison
#	     Ragothaman M. Yennamalli, Jaypee University of Information Technology
#Version: 1.0
#Last Modified: 26 Feb, 2019
#Description: The script filters out protein sequences with "X" amino acid.
#Input: Fasta filename containing protein sequences with unrecognized residues.
#Ouput: A fasta file with filtered out protein sequences
from Bio import SeqIO
import glob
import os
import sys
#
def genCleanFile(filename):
   path_name, file_name = os.path.split(filename)
   Sequences = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
   str2ryt_newFile = ""
   for key in Sequences:
      if ("X" not in Sequences[key].seq):
          str2ryt_newFile += Sequences[key].format("fasta")
   outfile = path_name+"/Cleaned_"+file_name
   with open(outfile, "w+") as f:
       f.write(str2ryt_newFile)
   f.close()
   return outfile
#

