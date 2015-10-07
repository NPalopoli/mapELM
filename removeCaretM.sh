#!/bin/bash
# Remove ^M at end of lines in UNIX files

sed -e "s///" "$1" >"$1".tmp
mv "$1".tmp "$1"
