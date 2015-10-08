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
from parseJPred import readFasta, readJPred

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

try:
  primaryAcc = sys.argv[2]
except IndexError:
  print("Primary accession code not specified. Format: ./mapELM.py <elm_instances[.date].tsv> <primaryAcc>")
  exit()

def readELMinstances(infile):
#  with open(infile, 'rb') as csvfile:
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
#  for row in elm:
#    print row
  return elm
#    for row in reader:
#    if row['Primary_Acc'] == 'Q63767':
#      print row['Primary_Acc'],row['Start'],row['End']

# Make dict of input files
parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
#predictions = readJPred(injnet)
#injnet.close()
#results = seq.copy()
#results.update(predictions)

def mapELMpositions(parsedELM,primaryAcc):
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = (row['Start'],row['End'])
#      ELMpos['Start'] = row['Start']
#      ELMpos['End'] = row['End']
  return ELMpos

ELMpos = mapELMpositions(parsedELM,primaryAcc)
print ELMpos

# Print table with results
#for row in zip(*([key] + value for key, value in results.items())):
#  print '\t'.join(map(str, row))

