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

# Setup

import os.path
import sys
from collections import OrderedDict
import csv
from Bio import SeqIO

# Read input files
try:
  inELMinstances = open(sys.argv[1])
#  inELMPDBFinddir = open(sys.argv[2])
  inELMPDBFinddir = sys.argv[2]
except IndexError:
  print("Input file(s) not specified. Format: ./getELMInstancesSeqs.py <elm_instances[.date].tsv> <ELMPDBFind_files>")  #<ELMPDBFind_parseSIFTS_<ELMID>.tsv>")
  exit()
except IOError:
  print("Input file(s) not found. Format: ./getELMInstancesSeqs.py <elm_instances[.date].tsv> <ELMPDBFind_parseSIFTS_<ELMID>.tsv>")
  exit()

# Functions

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

def readELMPDBFind(indir,infilename):
  '''Store ELM PDB mappings as list of dicts'''
  try:
    infile = open(os.path.sep.join([indir,infilename]))  # ELMPDBFind_parseSIFTS_ELMI002226.tsv
  except IOError:
    return
  parsedELMPDB = csv.DictReader(filter(lambda row: row[0]!='#' or row[0:13]=='#ELMAccession', infile),delimiter='\t', quotechar='"', )
  infile.close()
  return parsedELMPDB

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def printCSVDictReadertoFile(csvdictreader,outfile):
  '''Print list of dicts parsed with csv.DictReader to outfile'''
#  outfile = 'lala.outfile'
  with open(outfile,'wb') as fou:
    dw = csv.DictWriter(fou, delimiter='\t', fieldnames=csvdictreader.fieldnames, quoting=csv.QUOTE_ALL, lineterminator=os.linesep)
    headers = {}
    for n in dw.fieldnames:
      headers[n] = n
    dw.writerow(headers)
    for row in csvdictreader:
      dw.writerow(row)

def mapELMpositions(parsedELM,primaryAcc):
  '''Make dict with [ELMAccession:[Start,End,ELMType,ELMIdentifier,Primary_Acc]]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = [row['Start'],row['End'],row['ELMType'],row['ELMIdentifier'],row['Primary_Acc']]
  return ELMpos

def placeELM(seq,ELMpos,flankn,elmflankn):
  '''Map ELM to fasta sequence'''
  seqlen = len(seq['res'])
  seq['ELMpos'] = list('-' * seqlen)
  seq['ELMacc'] = list('-' * seqlen)
  seq['ELMType'] = list('-' * seqlen)
  seq['ELMIdentifier'] = list('-' * seqlen)
  seq['ELMflank'] = list('-' * seqlen)
#  seq['ELMflank10'] = list('-' * seqlen)
  flankn = flankn  # change number of flanking residues
  elmflankn = ''.join(['ELMflank',str(flankn)])
  seq[elmflankn] = list('-' * seqlen)
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
    # Compute flanking sequence
    flanksize = 20 + int(vals[0]) - 1 - int(vals[1])
    flanksizeodd = 0 if flanksize % 2 == 0 else 1
    flanksize = flanksize / 2
    flankfirst = max(1, int(vals[0]) - flanksize)
    flanklast = min(int(vals[1]) + flanksize + flanksizeodd, seqlen)
    if (flankfirst == 1 and flanklast == seqlen):
      continue
    elif (flankfirst == 1):
      flanklast = flanklast + flanksize - ( int(vals[0]) - 1 )
    elif (flanklast == seqlen):
      flankfirst = flankfirst - (int(vals[1]) + flanksize + flanksizeodd - seqlen)
    for pos in range(flankfirst-1,flanklast):
      seq['ELMflank'][pos] = seq['res'][pos]
    # Compute flanking 10 residues sequence
#    flankn = 10  # change number of flanking residues
    flanknfirst = max(1, int(vals[0]) - flankn)
    flanknlast = min(int(vals[1]) + flankn, seqlen)
    for pos in range(flanknfirst-1,flanknlast):
#      seq['ELMflank10'][pos] = seq['res'][pos]
      seq[elmflankn][pos] = seq['res'][pos]
  return seq

# Start

# Make dict of input files
#seq = readFasta(infasta)
#infasta.close()

parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
#ELMpos = mapELMpositions(parsedELM,primaryAcc)
#print parsedELM.writeheader()

outfile = 'elm_instances.2015-08-27_addPDB.tsv'
with open(outfile,'ab') as fou:
  dw = csv.DictWriter(fou, delimiter='\t', fieldnames=parsedELM.fieldnames, quoting=csv.QUOTE_ALL, lineterminator=os.linesep)
  headers = {}
  for n in dw.fieldnames:
    headers[n] = n
  dw.writerow(headers)
  for row in parsedELM:
    infilename = ''.join(['ELMPDBFind_parseSIFTS_',row['Accession'],'.tsv'])  # ELMPDBFind_parseSIFTS_ELMI002226.tsv
    parsedELMPDB = readELMPDBFind(inELMPDBFinddir,infilename)  # read ELMPDBFind_files_addheader/ELMPDBFind_parseSIFTS_ELMI002226.tsv
    if parsedELMPDB is not None:
#      print(parsedELMPDB)
      for pdb in parsedELMPDB:
#        print(pdb)
#        print(pdb.keys())
#        print(pdb.values())
        if pdb["PrimAccMatch"] == '-':
           continue
        if pdb["PrimAccMatch"] == 'Y':
#          row['PDB'] = row['PDB'] + ' ' + pdb["PDBID"]
          row['PDB'] = ' '.join([row['PDB'],pdb["PDBID"]])
    row['PDB'] = ' '.join(set(row['PDB'].split()))
    dw.writerow(row)
fou.close()
exit()

for row in parsedELM:
  infilename = ''.join(['ELMPDBFind_parseSIFTS_',row['Accession'],'.tsv'])  # ELMPDBFind_parseSIFTS_ELMI002226.tsv
  parsedELMPDB = readELMPDBFind(inELMPDBFinddir,infilename)  # read ELMPDBFind_files_addheader/ELMPDBFind_parseSIFTS_ELMI002226.tsv
  if parsedELMPDB is not None:
#    print(parsedELMPDB)
    for pdb in parsedELMPDB:
#      print(pdb)
#      print(pdb.keys())
#      print(pdb.values())
      if pdb["PrimAccMatch"] is 'Y':
        row['PDB'] = row['PDB'] + ' ' + pdb["PDBID"]
#      row['PDB'] = row['PDB'] + pdb["PDBID"]
#      print(row)
#  for pdb in parsedELMPDB:
#    print(pdb)
#    parsedELM['PDB'] = parsedELM['PDB'] + pdb['PDBID']
#    if (pdb['PrimAccMatch'] == 'Y'):
#      parsedELM['PDB'] = parsedELM['PDB'] + parsedELMPDB['PDBID']

#  outfilename = ''.join
#  printCSVDictReadertoFile(
#  inELMPDBFind.close()

#print(parsedELM)
#exit()

# dr.fieldnames contains values from first row of `f`.
outfile = 'elm_instances.2015-08-27_addPDB.tsv'
with open(outfile,'ab') as fou:
  dw = csv.DictWriter(fou, delimiter='\t', fieldnames=parsedELM.fieldnames, quoting=csv.QUOTE_ALL, lineterminator=os.linesep)
  headers = {}
  for n in dw.fieldnames:
    headers[n] = n
  dw.writerow(headers)
  dw.writerow(parsedELM)
  for row in parsedELM:
    dw.writerow(row)
fou.close()

exit()
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
'''
for instance in ELMpos:
  singleELMpos = {}
  singleELMpos[instance] = ELMpos[instance]
  flankn = 10  # set number of flanking residues
  elmflankn = ''.join(['ELMflank',str(flankn)])
  seq = placeELM(seq,singleELMpos,flankn,elmflankn)
  header = '|'.join([seq['name'][0],instance,ELMpos[instance][2],ELMpos[instance][3]])
  outname = instance + '.fasta'
  outfile=open(outname, 'w+')
  print >>outfile, '{}{}|sequence'.format('>',''.join(header))
  print >>outfile, ''.join(seq['res'])
  print >>outfile, '{}{}|instance'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMpos'])
  print >>outfile, '{}{}|flanks'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMflank'])
  print >>outfile, '{}{}|flank{}'.format('>',''.join(header),flankn)
  print >>outfile, ''.join(seq[elmflankn])
  outfile.close()
'''
