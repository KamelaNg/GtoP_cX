#!/usr/bin/env python
from os import system 
import glob 
import sys
import re
from collections import defaultdict
from sets import Set
from itertools import groupby
import argparse
import copy

"""
Author: Kamela Ng <13April2018>

This script converts whole genome sequence data in the form of tab files into the most liekly output that would be observed from a classicXpertMTB/RIF run on this sample

Input:
A folder containing all the tab files that are to be processed
A map file that states the relationship between codons where mutations are known to be detected by the classic XpertMTB/RIF (based on this following paper: Ng KC, Meehan CJ, Torrea G, Goeminne L, Diels M, Rigouts L, de Jong BC, Andre E. 2018. Potential Application of Digitally Linked Tuberculosis Diagnostics for Real-Time Surveillance of Drug-Resistant Tuberculosis Transmission: Validation and Analysis of Test Results. JMIR Med Inform 6:e12.) and positions in the genome
A filename for the output

The map file is optional. If not supplied, it is assumed the standard H37Rv NC000962.3 was used and that the file sample_mapfile.txt is in the current working folder
An example of this file is (tab separated columns):
Codon	Ref	Gen	Pos	
428	761088
430	761094
431	761097
432	761100
434	761106
435	761109
437	761115
441	761127
445	761139
446	761142
450	761154
452	761160

The filename for the output is optional. If not supplied, output will be genomeToProbe_Xpert.txt

Output:
The output file contains, for each input tab file:
Filename, RIF Resistance result, Mutant codon number, Mutant codon nucleotides, Absent classicXpertMTB/RIF probe, Observed reaction for each classicXpertMTB/RIF probe - 1 (present), 0 (absent)

Usage:
python classicXpert_GtoP_final.py --folder <folderName> --map <codonMapFile> (optional) --out <outputFilename> (optional)
"""
#Hardcode the codon to wildtype nucleotides
codon_nuc = {
			'428' : ['A', 'G', 'C'], 
			'430' : ['C', 'T', 'G'], 
			'431' : ['A', 'G', 'C'], 
			'432' : ['C', 'A', 'A'], 
			'434' : ['A', 'T', 'G'], 
			'435' : ['G', 'A', 'C'], 
			'437' : ['A', 'A', 'C'], 
			'441' : ['T', 'C', 'G'], 
			'445' : ['C', 'A', 'C'], 
			'446' : ['A', 'A', 'G'],
			'450' : ['T', 'C', 'G'], 
			'452' : ['C', 'T', 'G']
}

#hardcode the codon positions with the mutated codons and the associated probe change
classicXpertmut = {
			'428' : [['AGG'],'prA'],
			'430' : [['CCG'],'prA'], 
			'431' : [['GGC'],'prA'], 
			'432' : [['GAA'],'prB'],
			'434' : [['GTG','ATT','ACG'],'prB'], 
			'435' : [['GTC','TAC','TTC','GAA'],'prB'], 
			'437' : [['GAC'],'prC'], 
			'441' : [['CAG','TTG'],'prC'], 
			'445' : [['GGC','AGC','TCC','TAC','GAC','AAC','CAG','CAA'],'prD'], 
			'446' : [['CAG'],'prD'],
			'450' : [['TTG','TGG','TTC'],'prE'], 
			'452' : [['CCG'],'prE_del']}

#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True, help='Folder containing the tab files to be processed')
parser.add_argument('--map', required=False, default="sample_mapfile.txt", help='File listing the relationship between codon and genome position (default is "sample_maplfile.txt")')
parser.add_argument('--out', required=False, default="genomeToProbe_Xpert.txt", help='Filename for output (default is "genomeToProbe_Xpert.txt")')

args = parser.parse_args()
	
#read in map file 
try:
	mapf=open(args.map, 'rU')
except IOError:
	print "\n map file not found."
	sys.exit()

#Open output files
try:
	classicXpertOutputfinal = open(args.out, "w")  
except IOError:
	print 'no room for save file'						
	sys.exit()
classicXpertOutputfinal.write("Filename"+"\t"+"RIF Resistance"+"\t"+"MutCodonNum"+"\t"+"MutCodonVal"+"\t"+"AbsentProbe"+"\t"+"ProbeA"+"\t"+"ProbeB"+"\t"+"ProbeC"+"\t"+"ProbeD"+"\t"+"ProbeE"+"\t"+"ProbeE_del"+"\n")
  
#read in map and create a dictionary of the codons to genome positions
#each key is the three codon positions (supplied genome +0, +1, +2) and the value is the codon
mapdict={}
while 1:
	line=mapf.readline()
	if not line:
		break
	line=line.rstrip()
	if re.match("Codon",line):#on an info line so skip
		continue
	sections=line.split("\t")
	pos=sections[1]
	mapdict[pos, str(int(pos)+1), str(int(pos)+2)]=sections[0]
	codonpos=mapdict.keys() #get a list of all the genome positions, grouped as triplets
mapf.close()

#output probe pattern
outputList=["1","1","1","1","1","1"]
outputKeys=["prA","prB","prC","prD","prE","prE_del"]

#read in tab file
#for each file make a copy of the WT codon patterns, mutate it based on the vcf file and then output the associated result
OpenDir = glob.glob(args.folder + "/*")										  
for File in OpenDir: 
	#print File 
	if File.endswith(".tab"): 
		#print File
		try:
			tab=open(File,'rU')
			#print tab
		except IOError:
			print "\n tab file not found."
		sampleName=File.rsplit("_",8)[0].split("/")[-1]
		print 'Processing '+sampleName
		
		sampleCodons=copy.deepcopy(codon_nuc)
		while 1:
			line=tab.readline()
			#print s
			if not line:
				break
			line=line.rstrip()
			#print s
			if re.match("#",line):	#on an info line so skip
				continue
			sections=line.split("\t")
			#print sections
			mutpos=sections[0]
			altbase=sections[4]
			#print mutpos, altbase
			#go through the genome positions that are associated with any change in a codon base
			for sublist in codonpos:
				if mutpos in sublist: #check is the genome position from the vcf file is a position in a codon of interest
					codon=mapdict[sublist]
					#print codon
					if codon in classicXpertmut.keys(): #check if the codon is in the mutation list
						cXpertmutcodval = sampleCodons[codon]
						ind = sublist.index(mutpos)
						#modify the WT codon to be the new mutated codon
						cXpertmutcodval[ind] = altbase
						sampleCodons[codon]=cXpertmutcodval					
		tab.close()
			
		#output results
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		probeList=[]
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		codonList=codon_nuc.keys()
		codonList.sort()
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos]: #the codons are not the same 
				codonCheck="".join(sampleCodons[codonPos])
				if codonCheck in classicXpertmut[codonPos][0]: #if the mutated codon is a known one, get the probe change and change in the output
					probe=classicXpertmut[codonPos][1]
					sampleOutputList[outputKeys.index(probe)]="0"
					probeList.append(probe)
					mutCodonList.append(codonPos)
					mutCodonValList.append(codonCheck)
			
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			classicXpertOutputfinal.write(sampleName+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(probeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			classicXpertOutputfinal.write(sampleName+"\t"+"NOT DETECTED\n")
classicXpertOutputfinal.close()
sys.exit()