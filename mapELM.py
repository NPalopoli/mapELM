#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Parse elm_instances[.date].tsv from SLiMBench run
  File name: mapELM.py
  Author: Nicolas Palopoli
  Date created: 2015/10/05
  Date last modified: 2015/10/05
  Python Version: 2.7
'''

import sys
import csv
from collections import OrderedDict

# Read input files
try:
  inELMinstances = open(sys.argv[1])
#  injnet = open(sys.argv[2])
except IndexError:
  print("Input file not specified. Format: ./mapELM.py <elm_instances[.date].tsv>")
  exit()
except IOError:
  print("Input file not found. Format: ./mapELM.py <elm_instances[.date].tsv>")
  exit()

def readELMinstances(infile):
#  with open(infile, 'rb') as csvfile:
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  for row in elm:
    print row
  return elm
#    for row in reader:
#    if row['Primary_Acc'] == 'Q63767':
#      print row['Primary_Acc'],row['Start'],row['End']

"""
def readELMinstances(infile):
  '''Store ELM instances information from file.'''
  elm = OrderedDict()
#  readFirstSeq = False
  for line in infile:
#    if line[0] == '#' or line.strip() or line not in ['\n', '\r\n']:  # avoid header, empty or only whitespace lines
    if line[0] == '#':  # avoid header, empty or only whitespace lines
      continue
    line=line.rstrip()  # discard newline at the end (if any)
    words=line.split('\t')
    if words[0] == '"Accession"':  # store column names
      colNames = words  
    else:
      elm[words[0]] = words[1:]
#      name=words[0][1:]
#      seqs['res']=''
#    else :  # sequence, not header, possibly multi-line
#      seqs['res'] = seqs['res'] + line
#  seqs['res'] = list(seqs['res'])
#  seqs['position'] = range(1,len(seqs['res'])+1)
#  seqs['name'] = [name] * len(seqs['res'])
  print elm
  return elm
"""
# Make dict of input files
seq = readELMinstances(inELMinstances)
inELMinstances.close()
#predictions = readJPred(injnet)
#injnet.close()
#results = seq.copy()
#results.update(predictions)

print seq

# Print table with results
#for row in zip(*([key] + value for key, value in results.items())):
#  print '\t'.join(map(str, row))


'''
#ELM_Instance_Download_Version: 1.4
#ELM_Instance_Download_Date: 2015-08-28 00:05:52.031014
#Origin: elm.eu.org
#Type: tsv
#NumInstances: 2675
"Accession"	"ELMType"	"ELMIdentifier"	"ProteinName"	"Primary_Acc"	"Accessions"	"Start"	"End"	"References"	"Methods"	"InstanceLogic"	"PDB"	"Organism"
"ELMI002256"	"CLV"	"CLV_C14_Caspase3-7"	"ATN1_HUMAN"	"P54259"	"P54259 Q99495 Q99621 Q9UEK7"	"103"	"107"	"10085113 9535906"	"cleavage reaction; mutation analysis; western blot""true positive"	""	"Homo sapiens"
'''

"""
def readJPred(injnet):
  '''Store JPred predictions by program from file.'''
  entries={}  # dict for raw entries
  entries = OrderedDict()
  for line in injnet:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line = line.rstrip()  # discard newline at the end (if any)
      words = line.split(':')
      program = words[0]  # programs are keys
      entries[program] = ''
      values = words[1].split(',')  # predictions are values
      entries[program] = values[0:-1]
  tempValue = []
  for value in entries['JNETJURY']:  # replace empty values with dashes
    if value == '*':
      tempValue.append('*')
    else:
      tempValue.append('-')
  entries['JNETJURY'] = tempValue
  return entries
"""

'''
print seq
for name in seq:
  print name,len(seq[name]),list(seq[name])

for program in entries:
  print program,len(entries[program]),entries[program]
'''

'''
jnetpred
JNETCONF
JNETSOL25
JNETSOL5
JNETSOL0
JNETHMM
JNETPSSM
JNETJURY
JNETPROPE
JNETPROPH
JNETPROPC
'''

'''
# Split sequences by length
seqsSplit = {}  # dict for split seqs
for i in seqs.keys():
  chunk_size = 800  # max sequence length
  if len(seqs[i]) >= chunk_size:
    seqsSplit[i] = [seqs[i][pos:pos+chunk_size] for pos in range(0, len(seqs[i]), chunk_size)]
  else:
    seqsSplit[i] = list(seqs[i].split())

# Save fragments as separate files
for entry in seqsSplit.keys():
  count = 1
  for i in seqsSplit[entry]:
    outfile = entry+'_'+str(count)+'.fasta'  # output filename: 'fastaID_fragmentNumber.fasta'
    with open(outfile, 'a') as f:
      f.write('>{0}\n{1}\n'.format(entry,i))
    count += 1
'''
