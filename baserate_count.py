#!~/miniconda3/bin/python
from optparse import OptionParser


parser = OptionParser()
parser.add_option("--infile", dest="infile", help="give a .fa", metavar="FILE")
parser.add_option("--outfile",  dest="outfile", help="the name of oufput file [txt]", metavar="FILE")
parser.add_option('-i', dest="i", help='input')

(options, args) = parser.parse_args()

fqfile=open(options.infile,'r')
fqfilter=open(options.outfile,'w')


bases = ['A', 'C', 'G', 'T']


fqfilter.write("Base\t"+"A\t"+"T\t"+"C\t"+"G"+"\n")
for idx in range(0,int(options.i)):
	i=0
	base_dict = {'A':0,'T':0,'C':0,'G':0,'N':0}
	for line in open(options.infile,'r'):
		if i%2==1:
			base=line.strip()[idx].upper()
			base_dict[base] += 1
		i+=1		
	fqfilter.write(str(idx+1)+"\t"+str(base_dict["A"])+"\t"+str(base_dict["T"])+"\t"+str(base_dict["C"])+"\t"+str(base_dict["G"])+"\n")
		
fqfile.close()
fqfilter.close()
