from __future__ import division
import sys
import os

def calmp(num_of_reads, total_reads):
        mp = float(num_of_reads)/(int(total_reads)/1000000)
        return mp


def gettotalreadsfromflagstat(flagstatfile):
	f = open(flagstatfile)
        lines = f.readlines()
        mapped_reads = lines[0]
#        mapped_reads = int(x[0])
	print "mapped_reads", mapped_reads
	return mapped_reads


def main(Bedgraphfile, flagstatfile, readcountcorrBedgraphfile):
	total_reads = gettotalreadsfromflagstat(flagstatfile)
	wf = open(readcountcorrBedgraphfile, "w")
	f = open(Bedgraphfile)
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		num_of_reads = line[-1]
		mp = calmp(num_of_reads, total_reads)
		line[-1] = mp
		wf.write("\t".join(map(str,line))+"\n")
	wf.close()
	f.close()


if __name__=="__main__":
	Bedgraphfile, flagstatfile, readcountcorrBedgraphfile = sys.argv[1],sys.argv[2], sys.argv[3]
	main(Bedgraphfile, flagstatfile, readcountcorrBedgraphfile)
