# This script is like tomtom_finder.py, but for secondary structure motifs. It filters out RNA binding proteins that bind to the reverse
# complement of the input secondary structure motifs.

# eg. of use: python tomtom_finder_secondary.py  <TOMTOM output txt file> <TOMTOM output xml file>
# <output file to save info to>

#import necessary modules
import os
import sys

# import xml module
import xml.etree.ElementTree as ET

# read in input
input_tom_file = sys.argv[1]    #reads in input file containing results of TOMTOM operation (txt file)         
input_xml_file = sys.argv[2]    # reads in input file containing results of TOMTOM in xml

output_file = sys.argv[3]       # reads in output filename

tree = ET.parse(input_xml_file)
root = tree.getroot()

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

			# filter out ones that bind to the reverse complement
			if (line[9] == "+") :
				m.write(line[0] + "\t" + line[7] + "\t" + line[8] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + 
					protein + "\t" + line[2] + "\t" + line[6] + "\n")
				
m.close()











		

