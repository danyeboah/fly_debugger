#Daniel Yeboah-Kordieh
# This script reads in a set of sequences
# For each sequence, it generates a new random sequence based on the
# nucleotide frequencies of the input sequence.

# code to run each set of sequnces in RNApromo is commented out....to use this part of the code 
# uncomment  'cmdstring' lines. This must be run in the Cshell in a linux environment

# to use python random_seq_generator.py <file containing sequences> <no. of files containing random sequences to generate> 
# 

# import necessary modules
import sys
import random
import os

from Bio.Seq import Seq 
from Bio import SeqIO


# cmdstring1 = 'setenv SEGAL LAB Debug 1'
# os.system(cmdstring1) 

# read in file containing sequences
sequence_to_analyze = sys.argv[1]

# read in number of files to produce
no_files = int(sys.argv[2])

# initialize number of sequences to generate
sequence_no = 0

# read in output filename
#output_file = sys.argv[5]

# declare sequence variable
my_seq = Seq("") 
for seq_record in SeqIO.parse(sequence_to_analyze, "fasta"):
	my_seq = my_seq + seq_record.seq
	sequence_length = len(seq_record.seq)
	sequence_no += 1


Acount_ratio = float(my_seq.count('A')) / float(len(my_seq))
Gcount_ratio = float(my_seq.count('G')) / float(len(my_seq))
Ccount_ratio = float(my_seq.count('C')) / float(len(my_seq))
Tcount_ratio = float(my_seq.count('T')) / float(len(my_seq))


for n in range(0, no_files):
	#generate new filename
	filename = output_file + str(n + 1) + '.fasta'

	f = open(filename, 'w')
	
	for x in range(0, sequence_no):
		# sequence to store string letters, so they can be shuffled
		sequence_string_list = []

		for y in range(0, sequence_length):
			# generate a random number between 0 and 1
			random_no = random.uniform(0,1)


			# based on random number generated, generate random sequence
			if (random_no <= Acount_ratio):
				sequence_string_list.append('A') 
			elif (random_no <= Ccount_ratio):
				sequence_string_list.append('C') 
			elif (random_no <= Gcount_ratio):
				sequence_string_list.append('G') 
			else :
				sequence_string_list.append('T')


		# shuffle sequences
		random.shuffle(sequence_string_list)

		#build string from content of list
		sequence_string = ''.join(sequence_string_list)

		# write to file
		f.write('> ' + 'sequence_' + str(x + 1) + '\n')
		f.write(sequence_string + '\n')
		f.write('\n')

	# close file handle
	f.close()

#	cmdstring2 = 'perl rnamotifs08_motif_finder.pl -positive_seq ' + filename + ' -p n -output_pref random' + str(n + 1) + ' -output_dir random' + str(n + 1)

#	os.system(cmdstring2)





	

			





