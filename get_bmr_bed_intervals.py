import sys
import os
import subprocess
import commands

bedout = open("mm10_bmr_ref_intervals.bed", "w")
bmrs = [x.strip("\n") for x in open(sys.argv[1], "r")]
for x in range(len(bmrs)):
	###e.g. Scgb1b24:chr7:33738799-33750103
	info = commands.getoutput('head -n 1 mm10_bmr_references_06_28_2016/' + bmrs[x] + '.bmr')
	info_split = info.split(":")
	print(x)
	print(info_split)
	if len(info) == 0:
		continue
	
	chrom = info_split[1]
	pos = info_split[2].split("-")
	start = pos[0]
	end = pos[1]
	bedout.write("\t".join([chrom, start, end, bmrs[x]]) + "\n")


