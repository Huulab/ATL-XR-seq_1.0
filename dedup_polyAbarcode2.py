#!~/miniconda3/bin/python
from optparse import OptionParser
import pysam

parser = OptionParser()
parser.add_option("--inbam", dest="inbam", help="give a sorted bam to me", metavar="FILE")
parser.add_option("--outbam",  dest="outbam", help="the name of oufput bam", metavar="FILE")
(options, args) = parser.parse_args()

bamfile = pysam.AlignmentFile(options.inbam,'rb')
outbam = pysam.AlignmentFile(options.outbam, "wb", template=bamfile)



nums = list(range(1,23))+['X']+['M']
chrom = []
for i in nums:
	chrom.append('chr'+str(i))

for Chr in chrom:
	barcodes_hall=set()
	for read in bamfile.fetch(Chr):
		start_pos=str(read.reference_start)
		qname=read.query_name.split(":")[-1]
		qlength=str(read.query_length)
		barcodes=qname+"_"+start_pos+"_"+qlength
		if barcodes not in barcodes_hall:
			barcodes_hall.add(barcodes)
			outbam.write(read)


bamfile.close()
outbam.close()
