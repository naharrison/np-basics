#!/bin/csh -f

awk '{i++ ; n=(i%100); print $0 >> ("fileLists/files809810_" n ".list")}' fileLists/files809810.list
