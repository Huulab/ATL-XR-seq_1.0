#!~/miniconda3/bin/python
##discard 3 A base and screen no dipyrimidine
from optparse import OptionParser
import os


parser = OptionParser()
parser.add_option("--infile", dest="infile", help="give a cuted discared ployA fastq file to me", metavar="FILE")
parser.add_option("--outfile",  dest="outfile", help="the name of oufput file [fastq]", metavar="FILE")
(options, args) = parser.parse_args()


fqfile=open(options.infile,'r')
fqfilter=open(options.outfile,'w')

dipyrimidine=set(("TT","TC","CC","CT"))

i=0

for line in fqfile:
        #print i,i%4,line
	if i%4==0:
		seqID=line.strip("\n")
	elif i%4==1:
		sequence=line.strip("\n")
		dinucl=set()
		if len(sequence)>14 and len(sequence)<35:
			for pos in range(len(sequence)-1): 
				dinucl.add(sequence[pos:pos+2])
			if not dipyrimidine & dinucl:
				state="tooshort"
			else:
				if sequence[-2:len(sequence)]=="AA":
					state="AA"
					sequence=sequence[0:-2]
				elif sequence[-1]=="A":
					state="A"
					sequence=sequence[0:-1]
				else:
					state="normal"
					sequence=sequence
		else:
			state="tooshort"
	elif i%4==2:
		mark=line.strip("\n")
	elif i%4==3:
		qual=line.strip("\n")
		if state!="tooshort":
			if state=="AA":
				qual=qual[0:-2]
			elif state=="A":
				qual=qual[0:-1]
			else:
				qual=qual
			fqfilter.write(seqID+"\n"+sequence+"\n"+mark+"\n"+qual+"\n")
	i+=1

fqfile.close()
fqfilter.close()

