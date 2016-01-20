#!/bin/bash

awk 'BEGIN {FS="\t"} {if ($12!="\"\"") print $0}' /home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27.tsv >/home/npalopoli/SLiMBench/ELMmap/elm_instances.2015-08-27_withPDB.tsv
