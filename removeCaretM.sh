#!/bin/bash
# Remove ^M at end of lines in UNIX files

sed -e "s/
mv "$1".tmp "$1"