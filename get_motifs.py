# Daniel Yeboah-Kordieh
# Purpose: This script enables one to filter out the motif output of MEME based on the Branch length score and 
# the e-value in MEME
# This program reads in a file containing motif information from meme i.e. text file that is the result of meme
# It also reads in a file containing the evolutionary tree for the organism in question
# It then uses the information contained in the motif file and the evolutionary tree 
# to select significant motifs.
# This program requires that there is biopython installed on your machine

# It can be run from the command line using python get_motifs.py <newick file containing tree>
 #<file containing motifs> <output filename> eg. get_motifs.py drosophila_tree.nwk preliminary_motifs.txt motifs.txt
 
# import necessary modules
from Bio import Phylo
from Bio import motifs

# import sys
import sys

# import python math module
import math

# read in name of input newick file and motif file
inputnewick_file = sys.argv[1]
inputmotifs = sys.argv[2]

# read in name of output file
outputfile = sys.argv[3]

# read in name of mRNA
gene = sys.argv[4]

# open file
handle = open(inputmotifs)

# read in motif information
results = motifs.parse(handle, "meme")

# close file after reading in motif info
handle.close()

# create list to store motifs
motif_branch_lengths_list = []

# read in tree
tree = Phylo.read(inputnewick_file, "newick")

# create and print cutoff value
cutoff = (60/ (math.log10(0.01) + 100))
print ("cutoff = " + str(cutoff))

# function to calculate branch length of Drosophila with motifs
def motif_branch_length(list, tree, pvalue_correction):
	
	# find common ancestor of all species in list
	ancestor_all = tree.common_ancestor(list)
	
	# find distance of ancestor_all to all the species
	# use p value corrector to discount contribution of to total
	# branch length of motifs that are not very similar to consensus
	branch_length_motifs = 0
	for x in range(0, len(list)):
		branch_length_motifs += tree.distance(ancestor_all, list[x]) * pvalue_correction[x]
												
	# obtain sub corrections from branch_length_correction function
	sub_corrections = branch_length_correction(ancestor_all,list)

	# correct for branch length using corrections in sub_corrections list
	for correction in sub_corrections:
		branch_length_motifs -= correction	
	
	# divide branch length by total_branch_length to get branch length score
	# subtract 1 because distance of Drosophila genus to other genus is included which is 1
	total_length = tree.total_branch_length() - 1
	branch_length_score = branch_length_motifs/ total_length
	
	return branch_length_score

def branch_length_correction(ancestor, list):
	#  this function enables us to obtain corrections in the calculation of the branch
	#	length score.
	# The branch length score is first obtained by calculating the distance from the
	# most common ancestor of the species in question to the species. In doing so, there
	# will be double counts, so this function produces a list of the double counted figures
		
	for x in range(0, len(ancestor.clades)):
		# create list to store iterables of find any
		list_in_clade = []
		
		# for each clade, iterate to see how many of the species are in that clade
		for species in list:
			species_in_clade = ancestor.clades[x].find_any(species)
			if species_in_clade != None:
				list_in_clade.append(species_in_clade)
	
		#correct for double counting by subtracting distance to next ancestor if
		# there are more than 1 species in that clade		
		sub_correction_index = len(list_in_clade) - 1

		# where there is more than one species in clade, use a recursive algorithm
		# to correct for all multiple countings. Add all corrections to the sub_correction_array
		if sub_correction_index > 0:
			distance_to_ancestor = ancestor.distance(ancestor.clades[x])
			sub_correction = sub_correction_index * distance_to_ancestor
	
			sub_correction_array.append(sub_correction)
		
			branch_length_correction(ancestor.clades[x], list)

	return sub_correction_array

with open(outputfile, 'w') as f:
	f.write("cutoff = " + str(cutoff) + "\n")
	f.write("index" + "\t" + "BLS_score" + "\t" + "scaled_evalue" + "\t" + "motif_score" + 
	"\t" + "degenerate_consensus" + "\t" + "consensus" +  "\t" + "melano_sequence" + "\t" + "melano_start" +"\t" "melano_stop" + "\n")
f.closed
			
# iterate through each motif set and calculate branch length score for each
for x in range(0,len(results)):
	
	#list to store sub_corrections
	sub_correction_array = []
	
	# list to store p-values for each instance of the motif
	p_value_calc = []

	#list to store correction done were done using p values 
	# to affect contribution of each species to branch length score
	p_value_corrector = []
		
	species_name_list = []
	
	# iterate through instances of the motif, obtain pvalues of each instance and use
	# bonferroni holmes correction to correct for multiple hypothesis testing
	for y in range(0, len(results[x].instances)):
		p_value_calc.append(float(results[x].instances[y].pvalue))
	
	# calculate the log of the min pvalue for each motif
	max_pval = math.log10(min(p_value_calc))
	min_pval = math.log10(0.01)

	# convert each p-value to a 0-1 range using min and max pvalue
	for n in range(0, len(p_value_calc)):
		corrector = (math.log10(p_value_calc[n]) - min_pval)/(max_pval - min_pval)
		p_value_corrector.append(corrector)

	# iterate through instances of the motif, and for those instances with significant 
	# p values, calculate branch length
	for z in range(0, len(results[x].instances)):
		
		#obtain a list of species that have a motif in this instance
		species_name_list.append(results[x].instances[z].sequence_name)
	
	#calculate branch length for motif and add to motif_branch_lengths_list if melanogaster
	# is in the list
	if ("melanogaster" in species_name_list):
		
		BLS_score = (motif_branch_length(species_name_list, tree,p_value_corrector)) * 100

		scaled_evalue = math.log10(results[x].evalue) + 100   
		
		motif_score = BLS_score / scaled_evalue
		
		#obtain index of melanogaster in list for future use in obtaining start
		# and stop nucleotide of  motif for melanogaster
		melano_index = species_name_list.index("melanogaster")
		if (motif_score >= cutoff):
			with open(outputfile, 'a') as f:
				f.write(str(x+1) + "\t" + str(BLS_score) + "\t" + str(scaled_evalue)
				 + "\t" + str(motif_score) + "\t" + str(results[x].degenerate_consensus) + "\t"
				 + str(results[x].consensus) + "\t" + str(results[x].instances[melano_index]) + "\t"
				  + str(results[x].instances[melano_index].start) + "\t"
				   + str(results[x].instances[melano_index].start + results[x].instances[melano_index].length)+ "\n")
			#print weblogo
			results[x].weblogo(gene + "/motif_logo" + str(x + 1) + ".png")
			

