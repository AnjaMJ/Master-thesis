import sys,os,pysam,gzip
import numpy as np
from optparse import OptionParser
import subprocess

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input GWAS+freq file", default=None, type="string")
parser.add_option("-l", "--blockfile", dest="blockfile", help="Block file with PPA for each LD segment (def None)", default=None, type="string")
parser.add_option("-p", "--minppa", dest="minppa", help="Minimum block PPA allowed", default=0, type="float")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", default=None, type="string")
(options,args) = parser.parse_args()


blockfile = gzip.open(options.blockfile,"r")
infile = pysam.Tabixfile(options.infile, mode='r')
outfile = open(options.outfile,"w")
minppa = float(options.minppa)

readheader = infile.header
for x in readheader:
    header = x
header = header.split("\t")
header = list(filter(lambda a: a != "REF" and a != "ALT" and a != "ANC" and a != "DER" and a!= "SE" and a != "PPA", header))
header = "\t".join(header)
print >>outfile, header

blockheader = blockfile.readline()

for line in blockfile:
    regfields = line.strip("\n").split(" ")
    #print regfields
    regchr = regfields[2]
    regstart = int(regfields[3])
    regend = int(regfields[4])
    regppa = float(regfields[9])
    #print regchr,regstart,regend
    CurrLogPPA = -50
    CurrPPA = None
    CurrBest = None
    #try:
    for elem in infile.fetch(regchr,regstart-1,regend):
        fields = elem.strip().split("\t")
        #print fields
        PPA = fields[9]
        LogPPA = np.log10(float(fields[9]))
        if LogPPA > CurrLogPPA:
            CurrBest = fields[0:3]+[fields[7]]+fields[10:]
            CurrLogPPA = LogPPA
            CurrPPA = PPA
    #except:
    #    continue

    if regppa > minppa and CurrBest != None:

        print >>outfile, "\t".join(CurrBest)

#subprocess.call(["bgzip","-f",outfilename])
