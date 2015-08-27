# This script searches the output of TOMTOM (rna comparison tool) and outputs target motifs whose query motifs were deemed significant
# It reads in 3 files : 1. the file containing significant motifs from meme 2. file containing tomtom output in txt format
# and 3. file containing tomtom output in xml format
# It searches through the tomtom txt output for target motifs whose query motifs were significant according to get_motifs.py
# It then also searches the xml file to obtain the name of the rna binding protein for each significant query motif.

# eg. of use: python tomtom_finder.py <file containing info about significant motifs eg. gcl/gcl_motifs.txt> <TOMTOM output txt file> <TOMTOM output xml file>
# <output file to save info to>


#import necessary modules
import os
import sys


# import xml module
import xml.etree.ElementTree as ET

# read in input
input_seq_file = sys.argv[1]    #reads in input file containing selected motif sequences from MEME
input_tom_file = sys.argv[2]    #reads in input file containing results of TOMTOM operation (txt file)         
input_xml_file = sys.argv[3]    # reads in input file containing results of TOMTOM in xml

output_file = sys.argv[4]       # reads in output filename



tree = ET.parse(input_xml_file)
root = tree.getroot()

# create list to store indices of selected consensus sequences from MEME
index_consensus_seq = []

# obtain  consensus sequences from meme file and store in list
with open(input_seq_file) as seqs:
	for line in seqs:
		if not (line.startswith('cutoff') or line.startswith('index')):   # this condtion is to avoid reading the header lines
			line = line.rstrip("\n")
			line = line.split("\t")

			# store indices of consensus sequences in a list
			index_consensus_seq.append(line[0])






# open handle to write to output file and write headers
m = open(output_file, "w")
m.write("Query_ID" + "\t" + "Query_consensus" + "\t" + "Target_consensus" + "\t" + "p-value" + "\t" + 
	"e-value" + "\t" +"q-value" +"\t" "protein_name" + "\t" + "offset" + "\t" + "overlap" + "\n")


# select target sequence IDs and other info from only sequences that were found above and have right orientation
with open(input_tom_file) as toms:
	for line in toms:
		if not(line.startswith('#')):
			line = line.rstrip("\n")
			line = line.split("\t")


			for motif in root[1].iter('motif'):
				if (line[1] in motif.attrib['name']):
					protein = motif.attrib['alt']


			# only obtain info where consensus sequences is in the list of consensus sequences above
			if ((line[0] in index_consensus_seq) and (line[9] == "+")) :
				m.write(line[0] + "\t" + line[7] + "\t" + line[8] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + 
					protein + "\t" + line[2] + "\t" + line[6] + "\n")
				



m.close()











		

