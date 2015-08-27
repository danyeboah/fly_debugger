# This script is designed to figure out the 3'UTR region
# of an mRNA based on the mRNA's extended gene region and
# its coding region. It predicts the length of the the 3'UTR
# using the length of the known 3'UTR of a D. melanogaster orthologue
# It takes as input:
# 1. The name of folder containing the fasta files
# 2. Output filename

# It can be run from the command line using python three_prime_detector.py 
# <folder containing input fasta files> <output fasta filename>

# import necessary modules
import sys
import os
from Bio import SeqIO
from Bio import AlignIO

# read in input fasta foldername, and output filename from command line
inputfolder = sys.argv[1]
outputfile = sys.argv[2]


# create file to output data into
f = open(outputfile,'w')

melano_filepath = os.path.join(inputfolder, 'melanogaster.fasta')
melano_utr = SeqIO.read(melano_filepath, 'fasta')




# get utr_length. we assume that utr length is fairly conserved across species
utr_length = len(melano_utr) + 50

# write melanogaster three prime to final output file
f.write('> melanogaster')
f.write("\n")
for sequence in melano_utr.seq:
	f.write(sequence)
f.write("\n")
f.write("\n")

# for every file in folder that ends with .fasta 
# and does not begin with melanogaster, run clustalw2 on that file
for file in os.listdir(inputfolder):
	if (file.endswith('.fasta') and file.startswith('melanogaster') == False):
		filepath = os.path.join(inputfolder, file)
		cmdstring = 'clustalw2 ' + filepath
		os.system(cmdstring)

# for every file in folder that ends with .aln, read in the alignment
for file in os.listdir(inputfolder):
	if file.endswith('.aln'):
		filepath = os.path.join(inputfolder,file)
		
		# read in alignment
		alignment = AlignIO.read(filepath, 'clustal')
		
		# store the sequences as lists, so we can use list operations on them
		list1 = list(alignment[0])
		list2 = list(alignment[1])
		
		# look for index where three prime utr starts
		# reverse the list and iterate until you get the first letter
		
		for x in reversed(list2):
			if x!= '-':
				
				
				# to get index of the first letter in reverse, .index(x) will give you
				# the first occurence of x. You need the last occurence. Hence the code 
				# below
				threeprime_index_start = len(list2) - list2[::-1].index(x)
				threeprime_index_end = threeprime_index_start + utr_length + 1
				break
		
		# three prime utr is chunk of sequence within the range calculated on the 
		# extended gene region

		threeprime = alignment[0,threeprime_index_start:threeprime_index_end]
		
		# write the sequence to a fasta file
		fasta_id = '> ' + os.path.splitext(file)[0]
		f.write(fasta_id)
		f.write("\n")
		
		for sequencess in threeprime.seq:
			f.write(sequencess)
		f.write("\n")
		f.write("\n")
		
		
		
		
	