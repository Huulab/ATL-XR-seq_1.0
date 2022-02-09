#!~/miniconda3/bin/python
from optparse import OptionParser
import pysam

parser = OptionParser()
parser.add_option("--inbam", dest="inbam", help="give a fasta file.gz to me", metavar="FILE")
parser.add_option("--outfile",  dest="outfile", help="the name of oufput file [fasta]", metavar="FILE")
(options, args) = parser.parse_args()

bamfile = pysam.AlignmentFile(options.inbam,'rb')
fqfilter=open(options.outfile,'w')

nums = list(range(1,23))+['X']+['M']
chrom = []
for i in nums:
	chrom.append('chr'+str(i))
	
	

dict={}

for n in range(0,150):
	dict[str(n)]=0
	

for Chr in chrom:
	for read in bamfile.fetch(Chr):
		Anum=read.query_name.split(":")[-1][0:-1]
		dict[Anum]+=1

fqfilter.write("A length\t"+"count\n")
for a in range(10,101):
	fqfilter.write(str(a)+"\t"+str(dict[str(a)])+"\n")



bamfile.close()
fqfilter.close()
