#!/bin/bash

# Include PDBs found by ELMPDBFind_parseSIFTS.py in a copy of elm_instances.2015-08-27.tsv named elm_instances.2015-08-27_addPDB.tsv

head -5 /home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27.tsv >/home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27_addPDB.tsv

./elm_instances.2015-08-27_addPDB.py /home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27.tsv /home/npalopoli/20150924_ELM-Struct/ELMPDBFind/ELMPDBFind_files_addheader 1>elm_instances.2015-08-27_addPDB.out 2>elm_instances.2015-08-27_addPDB.err

