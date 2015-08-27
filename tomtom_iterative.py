# Daniel Yeboah-Kordieh
# This script reads in the results of an initial tomtom query with long motifs.
# It selects shorter motifs that were identified to be hits in the initial query and
# reruns them in tomtom.
# This is necessary because Tomtom compares motif with each other, but my aim is to find binding sites
# for motifs within the bigger motif. So the first tomtom query identifies this shorter motif and this script
# takes that motif from D. melanogaster and reruns it in Tomtom...that way I can use the statistical figures in 
# Tomtom more reliably

# To use run from command line with
# python tomtom_iterative.py <file containing motifs selected by get_motifs.py>
#							 <file containing result of Tomtom> <name of output meme file>



# import necessary modules
import os
import sys
from Bio.Seq import Seq


# import file containing selected motifs
input_motif_file = sys.argv[1]

# import filtered tomtom file
input_tomtom_file = sys.argv[2]

# read in output filename
output_meme_file = sys.argv[3]

# create associative array to store initial long motifs and spliced short motifs
long_motifs = {}
short_motifs = {}

# obtain the melanogaster sequences from the file containing selected motifs
with open(input_motif_file) as seqs:
	for line in seqs:
		if not (line.startswith('cutoff') or line.startswith('index')):   # this condtion is to avoid reading the header lines
			line = line.rstrip("\n")
			line = line.split("\t")

			long_motifs[line[0]] = line[6]

# slice longer motifs into shorter motifs based on info from tomtom
with open(input_tomtom_file) as seqs:
	
	#create list of keys....these keys are the motif indices
	dict_keys = []
	for line in seqs:
		if not (line.startswith('Query')):
			line = line.rstrip("\n")
			line = line.split("\t")
			dict_keys.append(line[0])
	
	# use keys in list to create dictionary 
	short_motifs = short_motifs.fromkeys(dict_keys)

	# add values for the keys
with open(input_tomtom_file) as seqs:
	for line in seqs:
		
		if not (line.startswith('Query')):
			line = line.rstrip("\n")
			line = line.split("\t")
			
			long_seq = long_motifs[line[0]]
			
			begin = abs(int(line[7]))
			length = int(line[8])

			split_motif = long_seq[begin: (begin + length)]

			# if there are no values associated with key, add value
			if (short_motifs[line[0]] == None):
				short_motifs[line[0]] = [split_motif]
						
			# only add to dictionary if it's not already in there
			elif (split_motif not in short_motifs[line[0]]):
				# add short motif to the dictionary, with the key being the index of the motif			
				short_motifs.setdefault(line[0],[]).append(split_motif)


# create command to be run in terminal
cmdstring = "./iupac2meme "
for item in short_motifs:
	for motif in short_motifs[item]:
		cmdstring = cmdstring + motif + " -named motif" + str(item) + " "

cmdstring = cmdstring + "-pseudo 1 " + ">" + output_meme_file

os.system(cmdstring)




