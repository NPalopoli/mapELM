#!/bin/bash

while read line; do echo "$line" | cut -f 5; done<elm_instances.2015-08-27.tsv | tail -n +7 |  tr -d '"' | sort -u >elm_instances.2015-08-27_PrimaryAcc.lst
