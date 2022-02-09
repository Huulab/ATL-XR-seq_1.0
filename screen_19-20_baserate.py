#!~/miniconda3/bin/python
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--infile", dest="infile", help="give a fa", metavar="FILE")
parser.add_option("--outfile",  dest="outfile", help="the name of oufput file [txt]", metavar="FILE")
parser.add_option("--t",  dest="dit", help="TT or TC")
parser.add_option("--n",  dest="n", help="position you want to know")
(options, args) = parser.parse_args()



fqfilter=open(options.outfile,'w')
pos=int(options.n)

bases = ['A', 'C', 'G', 'T']

count_1=[]
count_2=[]
count_3=[]
count_4=[]

with open(options.infile,'r') as fqfile:
	i=0
	for line in fqfile:
		if i%2==1:
			seq=line.strip()
			if len(seq)==26 and seq[pos-1:pos+1].upper()==options.dit:
				sequence=seq[pos-2:pos+2]
				count_1.append(sequence[0].upper())
				count_2.append(sequence[1].upper())
				count_3.append(sequence[2].upper())
				count_4.append(sequence[3].upper())
		i+=1

base_dict={}

for idx in range(0,4):
	base_dict[str(idx+1)]={}


for base in bases:
	base_dict["1"][base]=count_1.count(base)
	base_dict["2"][base]=count_2.count(base)
	base_dict["3"][base]=count_3.count(base)
	base_dict["4"][base]=count_4.count(base)
	
fqfilter.write("Base\t"+"A\t"+"T\t"+"C\t"+"G"+"\n")
	
		
for idx in range(0,4):
	fqfilter.write(str(idx+1)+"\t"+str(base_dict[str(idx+1)]["A"])+"\t"+str(base_dict[str(idx+1)]["T"])+"\t"+str(base_dict[str(idx+1)]["C"])+"\t"+str(base_dict[str(idx+1)]["G"])+"\n")

fqfile.close()
fqfilter.close()
