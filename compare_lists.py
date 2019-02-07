#!/usr/bin/env python

import sys


file1data = set([x.replace('\n','') for x in open(sys.argv[1],'r').readlines()])
file2data = set([x.replace('\n','') for x in open(sys.argv[2],'r').readlines()])
filename1 = sys.argv[1].split('/')[-1].split('.')[0]
filename2 = sys.argv[2].split('/')[-1].split('.')[0]

diff = file1data.difference(file2data)
overlap = file1data.intersection(file2data)

print filename1,'NOT',filename2,len(diff)
print filename1,'AND',filename2,len(overlap)

    