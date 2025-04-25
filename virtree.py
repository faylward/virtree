#!/usr/bin/env python
import sys, os, re, shlex, subprocess, pandas, argparse
from collections import defaultdict

#hmmdb = "hmm/vog_05p.hmm"

# run HMMER3
def run_hmmer(input_file, cpus, evalue, database, redo):
	output_file = re.sub(".faa", ".hmmout", input_file)
	cmd = "hmmsearch -E "+ str(evalue) +" --cpu "+ cpus +" --tblout "+ output_file +" "+ database +" "+ input_file	
	#print(cmd)
	cmd2 = shlex.split(cmd)
	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
		os.remove("out.txt")
	return output_file

# define function for parsing HMMER3 output
def parse_hmmout(hmmout):
	input = open(hmmout, "r")
	final_dict = defaultdict(int)
	hit_dict = defaultdict(lambda:"no_annot")
	bit_dict = defaultdict(float)

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			protein = tabs[0]
			hit = tabs[2]
			eval = float(tabs[4])
			score = float(tabs[5])
			if score > 30:
				if score > bit_dict[protein]:
					bit_dict[protein] = score
					hit_dict[protein] = hit
				else:
					pass
			else:
				pass
	for i in hit_dict:
		vog = hit_dict[i]
		final_dict[vog] +=1
	return final_dict

	

# main function that runs the program
def run_program(inputdir, project, database, minhit, evalue, cpus, redo):

	df = pandas.DataFrame()

	for i in os.listdir(inputdir):
		if i.endswith(".faa"):
			#name = re.sub("_genomic.fna.faa", "", i)
			inputfile = os.path.join(inputdir, i)
			#print(inputfile)
			#protein_file = predict_proteins(inputfile, inputdir)
			hmmout = run_hmmer(inputfile, cpus, evalue, database, redo)
			hit_dict = parse_hmmout(hmmout)
		#	fulllist = i.split(".")
		#	genomelist = fulllist[0:-2]
		#	genome = ".".join(genomelist)
			genome = re.sub(".faa", "", i)
			#print(genome)
			s1 = pandas.DataFrame(pandas.Series(hit_dict, name = genome))
			df = pandas.concat([df, s1], axis=1, sort=True)

	df = df.fillna(0)
	df2 = df.clip(0,1)

	rows_all_one = (df2==1).all(axis=1)
	#rows_to_drop = rows_all_one.tolist()
	#rows = [ind if j==1 in enumerate(rows_all_one.tolist()
	rows = [pair for pair in enumerate(rows_all_one.tolist()) if pair[1] == 1]
	rows_to_drop = [n[0] for n in rows]
	#print(rows_to_drop)
	#print(df2.shape)
	print("Removed the following VOGs because they were present in all genomes analyzed: ", list(df2.index[rows_to_drop]))
	if len(rows_to_drop) > 0:
		df2.drop(df2.index[rows_to_drop], inplace=True)
		#df2.drop(index=rows_to_drop, axis=1)
	else:
		pass

	#print(df2.shape)
	prof_tbl = project + "_profile.tsv"
	df2.to_csv(prof_tbl,sep="\t")
	#if BINARY == 1:
	#	df2 = df.clip(upper=1)
	#else:
	#	df2 = df.clip(upper=9)
	outputfile = project + "_bin.fna"
	o = open(outputfile, "w")
	for (columnName, columnData) in df2.items():
		vogsum = sum([float(n) for n in columnData])
		if vogsum >= minhit:
		#print(columnName, vogsum)
			string = "".join([str(int(n)) for n in columnData])
			o.write(">"+ columnName +"\n"+ string +"\n")



########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="hmmtree: Part of a workflow for making trees based on protein family presence/absence \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--inputdir', required=True, help='Input folder of protein FASTA files (ending in .faa)')
	args_parser.add_argument('-p', '--project', required=True, help='project name for outputs')
	args_parser.add_argument('-db', '--database', required=False, default="hmm/vog_05p.hmm", help='PATH to HMM database to use. Default is hmm/vog_05p.hmm (default), vog_025p.hmm, and vog_005p.hmm. See README for details')
	args_parser.add_argument('-g', '--minhit', required=False, default=int(1), help='minimum number of VOG hits that each viral region must have to be reported (default=4)')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-3), help='e-value that is passed to HMMER3 for the VOG hmmsearch (default=1e-3)')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser.add_argument('-v', '--version', action='version', version='ViralRecall v. 2.1')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	inputdir = args_parser.inputdir
	project = args_parser.project
	database = args_parser.database
	minhit = int(args_parser.minhit)
	evalue = str(args_parser.evalue)
	cpus = args_parser.cpus
	redo = args_parser.redo

	project = project.rstrip("/")

	run_program(inputdir, project, database, minhit, evalue, cpus, redo)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end




