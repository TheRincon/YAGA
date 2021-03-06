from Bio import Phylo
from ete2 import Tree
import argparse
from shutil import copy, copyfile
import os, errno
import sys
import json
from collections import Counter
import itertools
import ast
import random

# from Orthfinder actually get which are closer in distance to "target" and "trait"
def get_enriched_annotations(sco_list, r, og_file, target, directory, **others):
	final_array = []
	enriched_go = {}
	enriched_kegg = {}
	go_list = []
	total_go_list = []
	total_gene_list = []
	total_kegg_list = []
	# messy way to keep track of which header is target species, refine this
	# index = [i for i, s in enumerate(species_list) if x in s.lower()]
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
# i.e. so everyone can see what was used
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
		return c['ortho_finder_results_directory'], c['pop1'], c['pop2'], c["pop3"], c["pop4"]

def read_target_json_abba(tj):
	with open(tj, "rt") as tjf:
		j = tjf.read()
		c = json.loads(j)
		# protein_fasta_list = c.map(f: x, )
		return c['pop1_genome'], c['pop2_genome'], c["pop3_genome"], c["pop4_genome"], c["pop1_gff"], c["pop2_gff"], c["pop3_gff"], c["pop4_gff"]

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
# most likely interproscan doesn't have predictions for all genes
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
	return k

def cds_helper(genomes, gffs, directory):
	pop = [None] * 4
	dict_list = [None] * 4
	for i in range(0,4):
		uncombined = directory+"YAGA/GFFs/" + genomes[i].split("/")[-1].split(".")[0] + "_uncombined_cds.fasta"
		pop[i] = get_cds_coordinates(gffs[i], directory)
		run_bedtools(genomes[i], pop[i], uncombined)
		dict_list[i] = combine_cds_fastas(uncombined, directory)
	return dict_list

def run_bedtools(genome, gff, output):
	command_string = "bedtools getfasta -name -s -fi {} -bed {} > {}".format(genome, gff, output)
	os.system(command_string)

def combine_cds_fastas(fasta, directory):
	fasta_dict = {}
	current_id = ""
	with open(fasta, "rt") as ffasta:
		for line in ffasta:
			rows = line.split("\t")
			if current_id == rows[0][:-4] and rows[0][:-4] in fasta_dict:
				seq_line = next(ffasta)
				fasta_dict[rows[0][:-4]].append(seq_line.strip())
			else:
				#some cases the "getfasta" is opposite of what the -s (-) says
				old_id = current_id
				if not current_id == "":
					fasta_dict[current_id][-1] + "\n"
					if fasta_dict[current_id][-1].lower().startswith("atg") and not fasta_dict[current_id][0].lower().startswith("atg"):
						fasta_dict[current_id].reverse()
				current_id = rows[0][:-4]
				seq_line = next(ffasta)
				fasta_dict[rows[0][:-4]] = [seq_line.strip()]
	prefix = fasta.split("/")[-1].split(".")[0]
	pp = prefix.split("_")
	pop = "_".join(pp[:-2])
	with open(directory+"YAGA/GFFs/" + pop + "_final.fasta", "wt") as cdsfasta:
		for i, j in fasta_dict.iteritems():
			cdsfasta.write(i + "\n")
			for part in fasta_line_generator("".join(j).strip("\n")):
				cdsfasta.write(part + "\n")
	return (fasta_dict)

def orthogroup_mapping(species_map, species_list, orthogroups_file):
	pops_dict = {}
	orthogroup_map = []
	for x in species_map:
		pops_dict[x] = [i for i, s in enumerate(species_list) if x in s.lower()]
	with open(orthogroups_file, "rt") as fortho:
		next(fortho)
		for line in fortho:
			ogs = line.split("\t")
			line_dict = {}
			field_dict = {}
			if (len(ogs) > 4):
				for og in ogs[1:5]:
					for i, s in pops_dict.iteritems():
						field_dict[i] = ogs[int(s[0])+1]
				line_dict[ogs[0]] = field_dict
			orthogroup_map.append(line_dict)
	return (orthogroup_map, pops_dict)

def get_combinations(og_map, species_list):
	og_list_of_lists = {}
	final_combos = {}
	final = {}
	ucounter = 0
	lcounter = 0
	ulen = len(og_map)
	for g in og_map:
		ucounter = ucounter + 1
		sys.stdout.write('\r')
		sys.stdout.write("%d / %d" % (ucounter, ulen))
		sys.stdout.flush()
		x = ast.literal_eval(str(g))
		assert type(x) is dict
		for i, s in x.iteritems():
			combo = [[] for x in xrange(4)]
			y = ast.literal_eval(str(s))
			for k, v in y.iteritems():
				for h, m in species_list.iteritems():
					if (h == k):
						combo[int(m[0])] = list(v.split(","))
					og_list_of_lists[i] = combo
	print "\nCreating combinations..."
	for j, z in og_list_of_lists.iteritems():
		lcounter = lcounter + 1
		sys.stdout.write('\r')
		sys.stdout.write("%d / %d" % (lcounter, ulen))
		sys.stdout.flush()
		q = list(itertools.product(*z))
		final_combos[j] = q
	for a, b in final_combos.iteritems():
		fin = list(filter(lambda t: '' not in t, b))
		final[a] = fin
	return final

def get_seqs_for_alignments(combos, pop_files, combo_directory):
	alignment_prep = {}
	fcounter = 0
	no = len(combos)
	print "\nWriting combo fastas..."
	for x, s in combos.iteritems():
		fcounter = fcounter + 1
		sys.stdout.write('\r')
		sys.stdout.write("%d / %d" % (fcounter, no))
		sys.stdout.flush()
		ind = 0
		if len(s) > 10:
			truncated = s[:10]
		else:
			truncated = s
		for j in truncated:
			if (j[3] == '\r\n' or j[1] == '\r\n' or j[2] == '\r\n' or j[0] == '\r\n'):
				continue
			ind += 1
			with open(combo_directory + x + "_" + str(ind) + ".fasta", "wt") as fout:
				for y in pop_files:
					for q, r in y.iteritems():
						for z in range(0,4):
							if (">" + j[z].strip() == q):
								fout.write(q + "\n")
								pre_line = "".join(r).strip("\n")
								for parts in fasta_line_generator(pre_line):
									fout.write(parts + "\n")

def run_mafft(combo_directory, output_directory, ):
	mafft_list = os.listdir(combo_directory)
	xlen = len(mafft_list)
	xcounter = 0
	for mafft_file in mafft_list:
		xcounter + 1
		command_string = "mafft --quiet {} > {}".format(combo_directory + mafft_file, output_directory + "aligned_" + mafft_file)
		os.system(command_string)
		xcounter= xcounter + 1
		sys.stdout.write('\r')
		sys.stdout.write("%d / %d" % (xcounter, xlen))
		sys.stdout.flush()
	return len(mafft_list)

def fasta_line_generator(s):
	for start in range(0, len(s), 70):
		yield s[start:start+70]

def get_cds_coordinates(gff, directory):
	cds_list = []
	with open(gff, "rt") as gt:
		for line in gt:
			if line.startswith("#"):
				continue
			else:
				if line.split("\t")[2] == "CDS":
					parts = line.split("\t")
					descriptions = parts[8].split(";")
					for i in descriptions:
						if i.find("protein_id=") > -1:
							parts[2] = i[11:].strip()
							cds_list.append("\t".join(parts))
						elif i.find("transcript_id ") > -1:
							parts[2] = i[15:-1].strip()
							cds_list.append("\t".join(parts))
	with open(directory+"YAGA/GFFs/" + gff.split("/")[-1].split(".")[0] + "_cds.gff", "wt") as fout:
		for i in cds_list:
			fout.write(i)
	return directory+"YAGA/GFFs/" + gff.split("/")[-1].split(".")[0] + "_cds.gff"

def parse_abba_baba(abba_file, num, new_file, coarse_file):
	og_fine_dict = {}
	og_coarse_dict = {}
	with open(abba_file, "rt") as f:
		for line in f:
			clipped = line.strip().split("] ")
			if clipped[1].find("|||") > -1:
				og = next(f).strip().split("] ")[1][1:-1]
				subscript = next(f).strip().split("] ")[1][1:-1]
				next(f)
				next(f)
				next(f)
				D_1 = next(f).split("D statistic = ")[1]
				D = "".join(D_1[:-2])
				if D == 'NaN':
					D = "0.000"
				spacer = "_"
				og_fine_dict[og + "_" + subscript] = float(D)
				if og in og_coarse_dict:
					og_coarse_dict[og].append(float(D))
				else:
					og_coarse_dict[og] = [D]
	with open(new_file, "wt") as ffinal:
		for k, v in og_fine_dict.iteritems():
			ffinal.write("Orthogroup combintion: " + k + "\n")
			ffinal.write("D Statistic value = " + str(v) + "\n")
			ffinal.write("\n")

	#change D to 0.00 if negative?
	with open(coarse_file, "wt") as ffinal:
		for n, m in og_coarse_dict.iteritems():
			summy = 0.0
			for k in m:
				summy = summy + float(k)
			ffinal.write("Orthogroup: " + n + "\n")
			ffinal.write("D Statistic value = " + str(m) + "\n")
			ffinal.write("Average D value: " + str(summy / float(len(m))) + "\n")
			ffinal.write("\n")

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
		species_path = directory + "Orthogroups.csv"
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