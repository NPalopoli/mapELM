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
except IndexError:
  print("Input file not specified. Format: ./mapELM.py <elm_instances[.date].tsv> <primaryAcc>")
  exit()
except IOError:
  print("Input file not found. Format: ./mapELM.py <elm_instances[.date].tsv> <primaryAcc>")
  exit()

try:
  primaryAcc = sys.argv[2]
except IndexError:
  print("Primary accession code not specified. Format: ./mapELM.py <elm_instances[.date].tsv> <primaryAcc>")
  exit()

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def mapELMpositions(parsedELM,primaryAcc):
  '''Make dict with [ELMAccession:(Start,End)]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = (row['Start'],row['End'])
  return ELMpos

# Make list of dicts from input files
parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
# Print dict of start and end positions of ELMs for Primary Accession
ELMpos = mapELMpositions(parsedELM,primaryAcc)
print ELMpos
