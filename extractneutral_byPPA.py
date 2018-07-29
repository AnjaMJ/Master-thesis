import sys,os,pysam
from optparse import OptionParser
import subprocess

parser = OptionParser("$prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input GWAS+freq file", default=None, type="string")
parser.add_option("-p", "--maxppa", dest="maxppa", help="Maximum PPA allowed", default=1, type="float")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", default=None, type="string")
parser.add_option("-s", "--sep", dest="sep", help="Number of valid SNPs separating each printed SNP", default=None, type="int")
(options,args) = parser.parse_args()


infile = pysam.Tabixfile(options.infile, mode='r')
outfile = open(options.outfile,"w")

readheader = infile.header
for x in readheader:
    header = x
header = header.split("\t")
header = list(filter(lambda a: a != "REF" and a != "ALT" and a != "ANC" and a != "DER" and a!= "SE" and a != "PPA", header))
header = "\t".join(header)
print >>outfile, header

i=0
for line in infile.fetch():
    fields = line.strip().split("\t")
    PPA = float(fields[9])
    if PPA < options.maxppa:
        if i == options.sep:
            i = 0
            finalvec = fields[0:3]+[fields[7]]+fields[10:]
            print >>outfile, "\t".join(finalvec)
        i += 1
        
#subprocess.call(["bgzip","-f",outfilename])
