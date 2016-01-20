#!/bin/bash

# list all PDBs, one line per ELM instance
while read line; do echo "$line" | cut -f 12; done<elm_instances.2015-08-27_withPDB.tsv | tail -n +7 |  tr -d '"' | sort -u >elm_instances.2015-08-27_PDB.lst

# list all PDBs, one line per PDB
tr ' ' '\n' <elm_instances.2015-08-27_PDB.lst >elm_instances.2015-08-27_PDB-split.lst

# list all PDBs, one line per PDB, no duplicates
sort -u elm_instances.2015-08-27_PDB-split.lst >elm_instances.2015-08-27_PDB-split-nodups.lst
