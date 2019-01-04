from Bio import Phylo
from ete2 import Tree
import argparse
from shutil import copy, copyfile
import os, errno
import sys
import json
from collections import Counter

# from Orthfinder actually get which are closer in distance to "target" and "trait"
def get_enriched_annotations(sco_list, r, og_file, target, directory, **others):
	# too many lists, I actually don't even need all these
	final_array = []
	enriched_go = {}
	enriched_kegg = {}
	go_list = []
	total_go_list = []
	total_gene_list = []
	total_kegg_list = []
	# messy way to keep track of which header is target species, refine this
	c_count = 0
	h = 0
	with open(og_file, "rt") as f2:
		species_line_order = next(f2).split("\t")[1:]
		for i in species_line_order:
			c_count += 1
			if i.find(target) > -1:
				h = c_count
				# shaking my head
		for line in f2:
			tabs = line.split("\t")
			id  = tabs[0]
			# can I not just match from semantic name?
			target_set = tabs[h]
			if tabs[0] in sco_list:
				t_genes = target_set.split(",")
				# super messy
				if len(t_genes):
					for i in t_genes:
						total_gene_list.append(i)
						if i in r:
							enriched_go[i] = r[i]
						if ('kegg' in others):
							if i in kegg:
								enriched_kegg[i] = kegg[i]
	# pick better variable names: r, kk, gg
	for nn in r:
		gg = r[nn]
		if len(gg):
			for kk in gg:
				total_go_list.append(kk)
	for i in enriched_go:
		x = enriched_go[i]
		if len(x):
			for j in x:
				go_list.append(j)
	cos = Counter(total_go_list)
	cop = Counter(go_list)
	final_array.append("GOTERM" + "\t" + "ENRICHED_SET_AND_GO" + "\t" + "ALL_GOTERMS_MINUS_ENRICHED" + "\t" + "ENRICHED_NOT_IN_GO" + "\t" + "NOT_IN_GO_NOT_IN_ENRICHED" + "\n")
	for ff in cos:
		if ff in cop:
			final_array.append(str(ff) + "\t" + str(cop[ff]) + "\t" + str(int(cos[ff]) - int(cop[ff])) + "\t" + str(len(total_gene_list) - int(cop[ff])) + "\t" + str(10994 - (int(cos[ff]) - int(cop[ff])) - (len(total_gene_list) - int(cop[ff])) - int(cop[ff])) + "\n")
	del final_array[1]  # removes the '' row for GO terms, makes it look much nicer
	if ('kegg' in others):
		final_array.append("KEGGTERM" + "\t" + "ENRICHED_SET_AND_KEGG" + "\t" + "ALL_KEGGTERMS_MINUS_ENRICHED" + "\t" + "ENRICHED_NOT_IN_KEGG" + "\t" + "NOT_IN_KEGG_NOT_IN_ENRICHED" + "\n")
	with open(directory+"YAGA/ABBA_BABA_GO_OUTPUT.txt", "wt") as fi:
		for line in final_array:
			fi.write(line)

# 
def get_target_set(ts, lifestyle_match, outgroup, nearest_neighbor):
	tally = {}
	target_set = []
	for i, j in ts.items():
		if str(j).find(outgroup) > -1:
			tally[i] = outgroup
			# target_set.append(outgroup)
		elif str(j).find(nearest_neighbor) > -1:
			tally[i] = nearest_neighbor
			# target.append(nearest_neighbor)
		elif str(j).find(lifestyle_match) > -1:
			tally[i] = lifestyle_match
			target_set.append(i)
		else:
			tally[i] = "Tree appears to not be resolved"
			print "Error tree not resolved, Orthofinder failed to root this gene tree"
	return target_set


# extract GO terms from Interproscan Output
def get_go(txt):
	with open(txt, "rt") as f:
		current_id = ""
		go_dict = {}
		for line in f:
			rows = line.split("\t")
			transcript_id = rows[0]
			if current_id == rows[0] and len(rows) > 13 and rows[0] in go_dict:
				go_terms = rows[13].split("|")
				go_dict[rows[0]].extend(go_terms)
			elif current_id == rows[0] and len(rows) > 13:
				go_terms = rows[13].split("|")
				go_dict[rows[0]] = go_terms
			else:
				current_id = rows[0]
		for j in go_dict:
			# only get unique values for each list in the dictionary
			cc = unique(go_dict[j])
			go_dict[j] = cc
	return go_dict

def unique(list1):
	# I chose to make a new list here for readability, rather than a (lambda) filter
	unique_list = []
	for x in list1:
		if x not in unique_list:
			unique_list.append(x)
	return unique_list

# get KEGG from concatenated user_ko_text file if desired
def get_KEGG(g):
	kegg_dict = {}
	with open(g, "rt") as f3:
		for line in f3:
			parts = line.split("\t")
			gene = parts[0]
			if len(parts) > 1:
				ann = parts[1].strip()
				kegg_dict[gene] = ann
	return kegg_dict

# too sparse to contain useful info, no duplicates to look for enrichment
def get_cazy(v):
	cazy_dict = {}
	with open(v, "rt") as f4:
		for line in f4:
			parts = line.split("\t")
			gene = parts[0]
			if len(parts):
				ann = parts[1].strip()
				cazy_dict[gene] = ann
	return cazy_dict

# use ETE package to get leaves with "target" and get_sisters, build dict of sisters.
def get_neighbors(dst, species, target):
	num_of_species = len(species)
	target_sisters = {}
	headers = [[] for i in range(1, num_of_species)] # not even used, may use in the future
	for filename in os.listdir(dst):
		if filename.endswith(".txt"):
			with open(dst+filename, "rt") as f:
				content = f.read().strip()
				tree2 = Tree(content) # read the file in and make a tree. The way in the docs doesn't work for me
				for node in tree2.traverse("levelorder"): # could be any other way of traversing, just picked one
					if node.name.startswith(target): # barest semantic check possible
						x = node.get_sisters() # with MSA it never seems to return '', not so with "normal" Orthofinder
						target_sisters[filename.split("_")[0]] = node.get_sisters()
	return target_sisters

# safer call to sys exit then ending in main
def end_yaga():
	sys.exit()

# cannot distinguish between no space, already existing or errors, refine
def make_safe_dir(dir_name):
	try:
		os.makedirs(dir_name)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

# redundant, but duplicate single copy ortholog groups to a new directory for traceability
def copyfiles_to_dir(dst, src, sco_list):
	for og in sco_list:
		og_name = og.strip()
		copyfile(src + str(og_name) + "_tree.txt", dst + str(og_name) + "_tree.txt")

def find_names_species(file_species_overlaps):
	with open(file_species_overlaps, "rt") as f:
		species_line = f.readline()
		species = map(str.strip, species_line.split("\t")[1:])
		return species

# make a list to pass around scos, portability
def get_scos(tabs, orthogroups_file):
	og_array = []
	with open(orthogroups_file, "rt") as f1:
		next(f1)
		for line in f1:
			if line.find(tabs) > -1:
				og_array.append(line.split("\t")[0])
	return og_array

def get_missing():
	print "this"

def read_target_json(tj):
	with open(tj, "rt") as tjf:
		j = tjf.read()
		c = json.loads(j)
		# protein_fasta_list = c.map(f: x, )
		return c['target_species'], c['trait_match'], c["nearest_neighbor"], c["outgroup"]

# stupid but necessary, make sure DIRECTORY path ends with "/"
def append_slash(candidate_string):
	if type(candidate_string) != str:
		print "path must be a string"
		end_yaga()
	if not candidate_string.endswith("/"):
		candidate_string = candidate_string + "/"
		return candidate_string
	elif candidate_string.endswith("/"):
		return candidate_string

# not sure if this is Interproscan not annotating every gene or not, but I get 
# 10507 genes instead of the expected 10994 genes for Knufia. TODO look into this
def get_total_genes(txt):
	tot = []
	with open(txt, "rt") as f:
		c = next(f)
		v = c.split("\t")[0]
		tot.append(v)
		go_dict = {}
		for line in f:
			rows = line.split("\t")
			if v == rows[0]:
				continue
			else:
				tot.append(v)
				v = rows[0]
	k = unique(tot)

# return all the file paths I need in main from here, very messy
def directory_check(directory):
	end_path = ""
	middle_path = ""
	if directory.find("Results") < 0:
		print "Directory must be named like Results_MONTH[Day]_[COUNT] according to what Orthofinder generates"
		end_yaga()
	elif directory.find("Results") > -1:
		x = directory.find("Results")
		middle_path = directory[:x]
		end_path = directory[x+7:]
		end_path = append_slash(end_path)
		other_end = append_slash(directory[x+7:x+13])
		directory = append_slash(directory)
		species_path = directory + "Orthogroups_SpeciesOverlaps.csv"
		orthologues_path = directory + "Orthologues" + other_end
		return (middle_path, end_path, species_path, orthologues_path, directory, other_end)

# to make sure they are single copy I search the orthogroups file for ("\t1")^num + "\tnum"
# this is the expected pattern in the specification
def make_simple_tabs(num):
	tabs = ""
	for i in range(0, num):
		tabs = tabs + "\t1"
	tabs  = tabs + "\t" + str(num)
	return tabs