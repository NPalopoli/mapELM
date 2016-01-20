#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Get sequence of ELM Instances
  File name: getELMInstancesSeqs.py
  Author: Nicolas Palopoli
  Date created: 2016/01/19
  Date last modified: 2016/01/19
  Python Version: 2.7
'''

import sys
from collections import OrderedDict
import csv
from Bio import SeqIO

# Read input files
try:
  infasta = open(sys.argv[1])
  inELMinstances = open(sys.argv[2])
  primaryAcc = sys.argv[3]
except IndexError:
  print("Input file(s) not specified. Format: ./ELMotiSt.py <in.fasta> <elm_instances[.date].tsv> <primaryAcc>")
  exit()
except IOError:
  print("Input file(s) not found. Format: ./ELMotiSt.py <in.fasta> <elm_instances[.date].tsv> <primaryAcc>")
  exit()

def readFasta(infasta):
  '''Store fasta sequences from file.'''
  seqs = OrderedDict()
#  seqs={}  # dict for raw seqs
  readFirstSeq = False
  for line in infasta:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line=line.rstrip()  # discard newline at the end (if any)
      if line[0]=='>':  # or line.startswith('>'); distinguish header
        if readFirstSeq:  # exit if more than 2 sequences
          break
        readFirstSeq = True
        words=line.split()
        name=words[0][1:]
        seqs['res']=''
      else :  # sequence, not header, possibly multi-line
        seqs['res'] = seqs['res'] + line
  seqs['res'] = list(seqs['res'])
  seqs['position'] = range(1,len(seqs['res'])+1)
  seqs['name'] = [name] * len(seqs['res'])
  return seqs

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def mapELMpositions(parsedELM,primaryAcc):
  '''Make dict with [ELMAccession:[Start,End,ELMType,ELMIdentifier,Primary_Acc]]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = [row['Start'],row['End'],row['ELMType'],row['ELMIdentifier'],row['Primary_Acc']]
  return ELMpos

def placeELM(seq,ELMpos):
  '''Map ELM to fasta sequence'''
  seq['ELMpos'] = list('-' * len(seq['res']))
  seq['ELMacc'] = list('-' * len(seq['res']))
  seq['ELMType'] = list('-' * len(seq['res']))
  seq['ELMIdentifier'] = list('-' * len(seq['res']))
  seq['ELMflank'] = list('-' * len(seq['res']))
#  for accession, limits in ELMpos.iteritems():
#    for pos in range(int(limits[0])-1,int(limits[1])):
  for accession, vals in ELMpos.iteritems():
    for pos in range(int(vals[0])-1,int(vals[1])):
      seq['ELMpos'][pos] = seq['res'][pos]
      if '-' in seq['ELMacc'][pos]: 
        seq['ELMacc'][pos] = accession
        seq['ELMType'][pos] = vals[2]
        seq['ELMIdentifier'][pos] = vals[3]
      else:
        seq['ELMacc'][pos] = seq['ELMacc'][pos] + accession
        seq['ELMType'][pos] = seq['ELMType'][pos] + vals[2]
        seq['ELMIdentifier'][pos] = seq['ELMIdentifier'][pos] + vals[3]
    flanksize = ( 20 + int(vals[0]) - int(vals[1])) / 2
    flankfirst = max(0, int(vals[0]) - flanksize)
    flanklast = min(int(vals[1]) + flanksize, len(seq['res']))
    if (flankfirst == 0 and flanklast == len(seq['res'])):
      continue
    elif (flankfirst == 0):
      flanklast = flanklast + int(vals[0]) - 1
    elif (flanklast == len(seq['res'])):
#      flankfirst = flankfirst - (flanksize - (flanklast - len(seq['res'])))
      flankfirst = flankfirst - (flanklast - len(seq['res']))
    for pos in range(flankfirst-1,flanklast):
      seq['ELMflank'][pos] = seq['res'][pos]
  return seq

# Make dict of input files
seq = readFasta(infasta)
infasta.close()

parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
ELMpos = mapELMpositions(parsedELM,primaryAcc)
'''
print ELMpos
seq = placeELM(seq,ELMpos)

for instance in ELMpos:
  header = '|'.join([seq['name'][0],instance,ELMpos[instance][2],ELMpos[instance][3]])
print '{}{}|sequence'.format('>',''.join(header))
print ''.join(seq['res'])
#print '{}{}|{}|{}|{}|ELMInstance'.format('>',''.join(seq['name'][0]),ELMpos.keys(),ELMpos['ELMType'],ELMpos['ELMIdentifier'])
print ''.join(seq['ELMpos'])
'''
for instance in ELMpos:
  singleELMpos = {}
  singleELMpos[instance] = ELMpos[instance]
  seq = placeELM(seq,singleELMpos)
  header = '|'.join([seq['name'][0],instance,ELMpos[instance][2],ELMpos[instance][3]])
  outname = instance + '.fasta'
  outfile=open(outname, 'w+')
  print >>outfile, '{}{}|sequence'.format('>',''.join(header))
  print >>outfile, ''.join(seq['res'])
  print >>outfile, '{}{}|instance'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMpos'])
  print >>outfile, '{}{}|flanks'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMflank'])
  outfile.close()
