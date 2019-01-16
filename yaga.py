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
   abba     Gene by gene analysis for introgression (optionally MAUVE CoLinear Blocks)
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
		# prefixing the argument with -- means it's optional
		parser.add_argument('-d','--dir', help='Path to the Directory used to generate Orthofinder Results', required=True)
		parser.add_argument('-g','--go', help='Path to the a GO file (Interproscan results)', required=True)
		parser.add_argument('-k','--kegg', help='Path to the Kegg File (BLASTKoala file)', required=False)
		parser.add_argument('-t','--target', help='Path to the target json')
		args = vars(parser.parse_args(sys.argv[2:]))
		directory = args["dir"]
		go_file = args["go"]
		kegg_file = args["kegg"]
		tj = args["target"]
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
		target, lm, nn, out = yutils.read_target_json(tj)
		print "Examining tree structures and finding features..."
		ts = yutils.get_neighbors(directory+"YAGA/Single_Copy_Gene_Trees/", species, target)
		target_set = yutils.get_target_set(ts, lm, out, nn)
		c = yutils.get_total_genes(go_file)
		print "Writing output file..."
		yutils.get_enriched_annotations(target_set, gos, directory+"Orthogroups.csv", target, directory)
		command_string = "Rscript {} {} > {}".format(sys.path[0] + "/adjpvalue.r", directory + "YAGA/ABBA_BABA_GO_OUTPUT.txt", directory + "YAGA/FINAL_YAGA_OUTPUT.txt", directory + "YAGA/FINAL_YAGA_OUTPUT.txt")
		os.system(command_string)
		print "Finished!\nOutput in {}".format(directory+"YAGA/")

	# rename variables to pop1, pop2, etc. as these are more accepted terms
	def abba(self):
		print "Starting YAGA with \"abba\" option..."
		parser = argparse.ArgumentParser(description='Gene by gene analysis for introgression')
		parser.add_argument('-d','--dir', help='Path to the Directory used to generate Orthofinder Results', required=True)
		parser.add_argument('-t','--target', help='Path to the target json')
		parser.add_argument('-m', '--mauve', help='MAUVE alignment file', required=False)
		args = parser.parse_args(sys.argv[2:])
		directory = args.dir
		tj = args.target
		middle_path, end_path, species_path, orthologues_path, directory, other_end = yutils.directory_check(directory)
		print "Creating YAGA output directories...."
		yutils.make_safe_dir(directory+"YAGA/")
		yutils.make_safe_dir(directory+"YAGA/GFFs/")
		yutils.make_safe_dir(directory+"YAGA/Alignments/")
		yutils.make_safe_dir(directory+"YAGA/Combinations/")
		print "Getting Populations and Species..."
		species = yutils.find_names_species(species_path)
		pop1_genome, pop2_genome, pop3_genome, pop4_genome, pop1_gff, pop2_gff, pop3_gff, pop4_gff = yutils.read_target_json_abba(tj)
		print "Extracting CDS regions..."
		pop_dicts = yutils.cds_helper([pop1_genome, pop2_genome, pop3_genome, pop4_genome], [pop1_gff, pop2_gff, pop3_gff, pop4_gff], directory)
		pop1, pop2, pop3, pop4 = yutils.read_target_json(tj)
		og_dictionary, indexed_species = yutils.orthogroup_mapping([pop1, pop2, pop3, pop4], species, directory+"Orthogroups.csv")
		print "Generating ABBA/BABA Fastas from Orthogroups..."
		combinations = yutils.get_combinations(og_dictionary, indexed_species)
		yutils.get_seqs_for_alignments(combinations, pop_dicts, directory+"YAGA/Combinations/")
		print "Finished!\nOutput in {}".format(directory+"YAGA/")

	def baba(self):
		parser = argparse.ArgumentParser(description='Sliding window admixture test')
		parser.add_argument('--amend', action='store_true')
		args = parser.parse_args(sys.argv[2:])
		print 'Running abba, amend=%s' % args.amend

# still extremely messy TODO cleanup arg extraction, make more safety/sanity checks, look into
# protein2genome align, get DNA and get SNPs from there to do true ABBA/BABA analysis
if __name__ == "__main__":

	YAGA()