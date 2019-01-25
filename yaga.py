from Bio import Phylo
from ete2 import Tree
import argparse
from shutil import copy, copyfile
import os, errno
import sys
import json
from collections import Counter
import yutils

import gffutils
import pyfaidx

class YAGA(object):

	def __init__(self):
		parser = argparse.ArgumentParser(
			description='ABBA-BABA and Similar Analysis',
			usage='''yaga.py <command> [<args>]

The available commands are:
   yaga     Runs analysis on duplicated and missing orthogroups, and target species pairs
   abba     Gene by gene analysis for introgression
   baba     Sliding window admixture test

''')

		parser.add_argument('command', help='Subcommand to run')
		# parse_args defaults to [1:] for args, but you need to
		# exclude the rest of the args too, or validation will fail
		args = parser.parse_args(sys.argv[1:2])
		if not hasattr(self, args.command):
			print 'Unrecognized command'
			parser.print_help()
			exit(1)
		# use dispatch pattern to invoke method with same name
		getattr(self, args.command)()

	def yaga(self):
		print "Starting YAGA with \"yaga\" option..."
		parser = argparse.ArgumentParser(description='Record changes to the repository')
		# parser.add_argument('-d','--dir', help='Path to the Directory used to generate Orthofinder Results', required=True)
		parser.add_argument('-g','--go', help='Path to the a GO file (Interproscan results)', required=True)
		parser.add_argument('-k','--kegg', help='Path to the Kegg File (BLASTKoala file)', required=False)
		parser.add_argument('-t','--target', help='Path to the target json', required=True)
		args = vars(parser.parse_args(sys.argv[2:]))
		go_file = args["go"]
		kegg_file = args["kegg"]
		tj = args["target"]
		directory, pop1, pop2, pop3, pop4 = yutils.read_target_json(tj)
		middle_path, end_path, species_path, orthologues_path, directory, other_end = yutils.directory_check(directory)
		species = yutils.find_names_species(species_path)
		sco_list = yutils.get_scos(yutils.make_simple_tabs(len(species)), directory + "Orthogroups.GeneCount.csv")
		src = orthologues_path + "Gene_Trees/"
		gos = yutils.get_go(go_file)
		if args["kegg"] is not None:
			x = yutils.get_KEGG(kegg_file)
		print "Creating YAGA output diretories...."
		yutils.make_safe_dir(directory+"YAGA/")
		yutils.make_safe_dir(directory+"YAGA/Single_Copy_Gene_Trees/")
		print "Copying Single copy orthogroup files to \"/YAGA/Single_Copy_Gene_Trees/\"..."
		yutils.copyfiles_to_dir(directory+"YAGA/Single_Copy_Gene_Trees/", src, sco_list)
		print "Examining tree structures and finding features..."
		ts = yutils.get_neighbors(directory+"YAGA/Single_Copy_Gene_Trees/", species, target)
		target_set = yutils.get_target_set(ts, pop3, pop4, pop2)
		c = yutils.get_total_genes(go_file)
		print "Writing output file..."
		yutils.get_enriched_annotations(target_set, gos, directory+"Orthogroups.csv", target, directory)
		command_string = "Rscript {} {} > {}".format(sys.path[0] + "/adjpvalue.r", directory + "YAGA/ABBA_BABA_GO_OUTPUT.txt", directory + "YAGA/FINAL_YAGA_OUTPUT.txt", directory + "YAGA/FINAL_YAGA_OUTPUT.txt")
		os.system(command_string)
		print "Finished!\nOutput in {}".format(directory+"YAGA/")

	def abba(self):
		print "Starting YAGA with \"abba\" option..."
		parser = argparse.ArgumentParser(description='Gene by gene analysis for introgression')
		parser.add_argument('-t','--target', help='Path to the target json', required=True)
		# parser.add_argument('-m', '--mauve', help='MAUVE alignment file', required=False)
		args = parser.parse_args(sys.argv[2:])
		tj = args.target
		directory, pop1, pop2, pop3, pop4 = yutils.read_target_json(tj)
		middle_path, end_path, species_path, orthologues_path, directory, other_end = yutils.directory_check(str(directory))
		print "Creating YAGA output directories...."
		yutils.make_safe_dir(directory+"YAGA/")
		yutils.make_safe_dir(directory+"YAGA/GFFs/")
		yutils.make_safe_dir(directory+"YAGA/Alignments/")
		yutils.make_safe_dir(directory+"YAGA/Combinations/")
		yutils.make_safe_dir(directory+"YAGA/Results/")
		print "Getting Populations and Species..."
		species = yutils.find_names_species(species_path)
		pop1_genome, pop2_genome, pop3_genome, pop4_genome, pop1_gff, pop2_gff, pop3_gff, pop4_gff = yutils.read_target_json_abba(tj)
		print "Extracting CDS regions..."
		pop_dicts = yutils.cds_helper([pop1_genome, pop2_genome, pop3_genome, pop4_genome], [pop1_gff, pop2_gff, pop3_gff, pop4_gff], directory)
		og_dictionary, indexed_species = yutils.orthogroup_mapping([pop1, pop2, pop3, pop4], species, directory+"Orthogroups.csv")
		print "Generating ABBA/BABA Fastas from Orthogroups..."
		combinations = yutils.get_combinations(og_dictionary, indexed_species)
		yutils.get_seqs_for_alignments(combinations, pop_dicts, directory+"YAGA/Combinations/")
		print "\nGetting MAFFT alignments..."
		num_of_alignments = yutils.run_mafft(directory+"YAGA/Combinations/", directory+"YAGA/Alignments/")
		print "\nCalling into R to calculate ABBA/BABA likelihoods..."
		command_string = "Rscript {} {} > {}".format(sys.path[0] + "/abbababa.r", directory + "YAGA/Alignments/", directory + "YAGA/Results/ABBA_BABA_OUTPUT.txt")
		os.system(command_string)
		yutils.parse_abba_baba(directory+"YAGA/Results/ABBA_BABA_OUTPUT.txt", 9418, directory+"/YAGA/Results/Final_ABBA_BABA.txt", directory+"/YAGA/Results/Final_Orthogroup_ABBA_BABA.txt") # yutils.parse_abba_baba(directory+"YAGA/Results/ABBA_BABA_OUTPUT.txt", num_of_alignments)
		print "Finished!\nOutput in {}".format(directory+"YAGA/")

	def baba(self):
		print "Starting YAGA with \"baba\" option..."
		parser = argparse.ArgumentParser(description='Sliding window admixture test')
		parser.add_argument('-d','--dir', help='Path to the Directory used to generate Orthofinder Results', required=True)
		parser.add_argument('-m', '--mauve', help='MAUVE alignment file', required=False)
		parser.add_argument('-t','--target', help='Path to the target json', required=True)
		args = parser.parse_args(sys.argv[2:])
		tj = args.target
		directory = args.dir
		middle_path, end_path, species_path, orthologues_path, directory, other_end = yutils.directory_check(directory)

if __name__ == "__main__":

	YAGA()