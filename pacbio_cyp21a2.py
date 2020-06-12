#!/usr/bin/env python
# encoding: utf-8

""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       pacbio_cyp21a2.py                                                  ///
 ///                                                                          ///
///////     Structure detection for CYP21A2 / CYP21A1P rearrangements        //////                                                                 //////
   ///      using PacBio long reads                                            ///
  ///                                                                         ///
 ///        Written by:     Zach Stephens                                    ///
///////     For:            Health Sciences Research, Mayo Clinic, MN       ///////
   ///      Date:           July 26, 2018                                      ///
  ///       Contact:        zstephe2@illinois.edu                             ///
 ///                                                                         ///
/////////////////////////////////////////////////////////////////////////////// """

from __future__ import print_function

import os
import re
import sys
import copy
import gzip
import time
import argparse
import subprocess
import numpy as np
import pickle

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1])
sys.path.append(SIM_PATH+'/py/')

from pbcyp_graphFunc import exhaustive_DAG
from pbcyp_misc      import exists_and_is_nonZero, makedir, getInputFileExtension, RC, read_fq_entry, read_fa_entry, print_nicely, l_2_cd, lcs
from pbcyp_plotFunc  import plotting_object
from pbcyp_vcf       import writeVCF


"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""


parser = argparse.ArgumentParser(description='pacbio_cyp21a2.py')
parser.add_argument('-i',  type=str,   required=True,  metavar='<str>',                    help="* input.fq")
parser.add_argument('-o',  type=str,   required=True,  metavar='<str>',                    help="* output_dir/")
parser.add_argument('-m',  type=str,   required=False, metavar='<str>',     default='',    help="motifs.p")
parser.add_argument('-l',  type=int,   required=False, metavar='<int>',     default=500,   help="minimum read length")
parser.add_argument('-p',  type=str,   required=False, metavar='<str>',     default='',    help="preprocessed_reads.p")
# "nice motif pairs" parameters
parser.add_argument('-ng', type=int,   required=False, metavar='<int>',     default=1,     help="(nice pairs) max contiguous motif gap")
parser.add_argument('-nd', type=int,   required=False, metavar='<int>',     default=5,     help="(nice pairs) max distance for motif gaps")
parser.add_argument('-np', type=float, required=False, metavar='<float>',   default=0.30,  help="(nice pairs) palindrome threshold")
# quality score thresholds
parser.add_argument('-qm', type=int,   required=False, metavar='<int>',     default=15,    help="minimum quality score for motif")
parser.add_argument('-qp', type=int,   required=False, metavar='<int>',     default=30,    help="minimum quality score for motif (prio)")
parser.add_argument('-qc', type=float, required=False, metavar='<int>',     default=0.10,  help="minimum percentage of read length explained by motifs")
# graph scoring parameters
parser.add_argument('-gf', type=int,   required=False, metavar='<int>',     default=10,    help="max gap for free connections")
parser.add_argument('-go', type=int,   required=False, metavar='<int>',     default=200,   help="permitted overlap (bp) in motif spans")
parser.add_argument('-gm', type=int,   required=False, metavar='<int>',     default=1000,  help="max length of unexplained read between motifs")
# clustering parameters
parser.add_argument('-cc', type=int,   required=False, metavar='<int>',     default=100,   help="max gap for joining distant motifs")
parser.add_argument('-cd', type=int,   required=False, metavar='<int>',     default=10,    help="max motif distance for clustering")
parser.add_argument('-cf', type=float, required=False, metavar='<float>',   default=0.50,  help="candidate structure must be close to at least this fraction of reads in an existing cluster to join it")
# plotting parameters
parser.add_argument('-pr', type=int,   required=False, metavar='<int>',     default=0,     help="plot this many raw read data (for debugging)")
parser.add_argument('-pc', type=int,   required=False, metavar='<int>',     default=100,   help="plot clusters supported by at least this many reads")
parser.add_argument('-pa', type=int,   required=False, metavar='<int>',     default=10,    help="max #motifs from edge to be considered anchored (high conf)")
parser.add_argument('--extra-plots',   required=False, action='store_true', default=False, help='produce extra plots (histograms, weak clusters)')
parser.add_argument('--skip-qs',       required=False, action='store_true', default=False, help='skip read quality score filter')
args = parser.parse_args()

# basic parameters
(INF, OUT_DIR, IN_MOTIFS, OUF_READS, MIN_RLEN) = (args.i, args.o, args.m, args.p, args.l)

if OUT_DIR[-1] == '/':
	OUT_DIR += '/'
makedir(OUT_DIR)

# read motif data
if IN_MOTIFS == '':
	print('using default motifs: cyp21a2-a1p-motifs.p')
	IN_MOTIFS = SIM_PATH+'/cyp21a2-a1p-motifs.p'
[MOTIFS, MOTIF_ALTS, N_PRIMER_MOTIFS, MOTIF_TRIM, PLOTTING_META] = pickle.load(open(IN_MOTIFS,"rb"))

# "nice motif pairs" parameters
(GAP_FLEX, DIST_WIGGLE, PALINDROME_PERCENTAGE) = (args.ng, args.nd, args.np)

# hard-coded values for PacBio quality score scheme
PACBIO_QSCORE_OFFSET = 33
PACBIO_QSCORE_MAX    = 100
# quality score thresholds
(QSCORE_MIN_ALL_MOTIF, QSCORE_MIN_PRIO_POS, MIN_READ_COVERAGE) = (args.qm, args.qp, args.qc)
(SKIP_QUALSCORE_FILTERS) = (args.skip_qs)

# graph parameters
(GRAPH_MAX_FREE_GAP, ALLOWED_OVERLAP_MAX, ALLOWED_GAP_MAX) = (args.gf, args.go, args.gm)

# clustering parameters
(JOIN_CONTIGUOUS, MAX_SHIFT_DIST, MAX_SHIFT_DIST_FRAC) = (args.cc, args.cd, args.cf)

# plotting parameters
(N_READS_TO_PLOT, MIN_READS_PER_CLUSTER, HIGH_CONF_ANCHOR) = (args.pr, args.pc, args.pa + N_PRIMER_MOTIFS)
(EXTRA_PLOTS) = (args.extra_plots)
my_plotting = plotting_object(MOTIFS,PLOTTING_META)

# parameters used for deriving motifs (to be included in output.p, if specified)
PARAMETERS = {'GAP_FLEX':GAP_FLEX,
              'MIN_RLEN':MIN_RLEN,
              'MOTIF_TRIM':MOTIF_TRIM,
              'DIST_WIGGLE':DIST_WIGGLE,
              'PACBIO_QSCORE_OFFSET':PACBIO_QSCORE_OFFSET,
              'PACBIO_QSCORE_MAX':PACBIO_QSCORE_MAX,
              'QSCORE_MIN_ALL_MOTIF':QSCORE_MIN_ALL_MOTIF,
              'QSCORE_MIN_PRIO_POS':QSCORE_MIN_PRIO_POS}

#
FUZZY = 1
if FUZZY:
	import regex

# motif character assignment
#
ASCII_OFFSET = 33
DEAD_CHAR    = chr(ASCII_OFFSET-2)
symbols = [[],[]]
decode  = {}
for i in range(len(MOTIFS)):
	symbols[0].append(chr(ASCII_OFFSET + i))
	decode[chr(ASCII_OFFSET + i)] = 'G_'+str(i)
for i in range(len(MOTIFS)):
	symbols[1].append(chr(ASCII_OFFSET + len(MOTIFS) + i))
	decode[chr(ASCII_OFFSET + len(MOTIFS) + i)] = 'P_'+str(i)

# construct all possible "nice motif pairs" that we'll be looking for in the reads
#
linked_motif_pairs = []
for i in range(len(MOTIFS)-1):
	if (i < N_PRIMER_MOTIFS and i+1 < N_PRIMER_MOTIFS) or (i >= len(MOTIFS)-N_PRIMER_MOTIFS and i+1 >= len(MOTIFS)-N_PRIMER_MOTIFS):
		continue
	linked_motif_pairs.append((symbols[0][i], symbols[0][i+1], MOTIFS[i+1][2] - MOTIFS[i][2], DEAD_CHAR, []))
	linked_motif_pairs.append((symbols[1][i], symbols[1][i+1], MOTIFS[i+1][3] - MOTIFS[i][3], DEAD_CHAR, []))
for gf in range(GAP_FLEX):
	for i in range(len(MOTIFS)-2-gf):
		if (i < N_PRIMER_MOTIFS and i+2+gf < N_PRIMER_MOTIFS) or (i >= len(MOTIFS)-N_PRIMER_MOTIFS and i+2+gf >= len(MOTIFS)-N_PRIMER_MOTIFS):
			continue
		linked_motif_pairs.append((symbols[0][i], symbols[0][i+2+gf], MOTIFS[i+2+gf][2] - MOTIFS[i][2], ''.join(symbols[0][i+1:i+2+gf]), [MOTIFS[i+1+n][2] - MOTIFS[i][2] for n in range(gf+1)]))
		linked_motif_pairs.append((symbols[1][i], symbols[1][i+2+gf], MOTIFS[i+2+gf][3] - MOTIFS[i][3], ''.join(symbols[1][i+1:i+2+gf]), [MOTIFS[i+1+n][3] - MOTIFS[i][3] for n in range(gf+1)]))


"""//////////////////////////////////////////////////
////////////          READ INPUT         ////////////
//////////////////////////////////////////////////"""


count_min_all  = {n:0 for n in range(PACBIO_QSCORE_MAX+1)}
count_min_pri  = {n:0 for n in range(PACBIO_QSCORE_MAX+1)}
discard_reads1 = 0
discard_reads2 = 0
discard_reads3 = 0
inParams       = {}

if not exists_and_is_nonZero(INF):
	print('\nError: Input file does not exist or is empty.\n')
	exit(1)

(isPICKLE, isFASTQ, isFASTA, isGZIP) = getInputFileExtension(INF)

if isPICKLE:
	tt = time.time()
	sys.stdout.write('reading input pickle... ')
	sys.stdout.flush()
	[inParams, all_data] = pickle.load(open(INF,"rb"))
	print(int(time.time()-tt),'(sec)')

elif isFASTQ or isFASTA:

	if isFASTQ:
		lineDivisor = 4
		message_out = 'reading input fastq'
	else:
		lineDivisor = 2
		message_out = 'reading input fasta'

	# get read count from wc -l if not gzipped, gzcat --> wc -l if gzipped...
	# #### is hardcoding gzcat ok? or would zcat be used instead?
	if not isGZIP:
		(out,err) = subprocess.Popen(['wc','-l',INF], stdout=subprocess.PIPE).communicate()
		nTOT      = int(int([n for n in out.decode('utf-8').split(' ') if len(n)][0])/lineDivisor)
		f = open(INF,'r')
	else:
		ps   = subprocess.Popen(['gzcat',INF], stdout=subprocess.PIPE)
		out  = subprocess.check_output(['wc','-l'], stdin=ps.stdout)
		ps.wait()
		nTOT = int(int([n for n in out.decode('utf-8').split(' ') if len(n)][0])/lineDivisor)
		f = gzip.open(INF,'rt')
		message_out += ' (gzipped)'

	nReads   = 0
	tt       = time.time()
	all_data = []
	print(message_out + '...')

	# read all the reads...
	while True:
		if isFASTQ:
			(my_read, my_qual, my_rname) = read_fq_entry(f)
		else:
			(my_read, my_qual, my_rname) = read_fa_entry(f)
		if not len(my_read):
			break

		my_qual = [min([ord(n) - PACBIO_QSCORE_OFFSET, PACBIO_QSCORE_MAX]) for n in my_qual]
		rdat_it = [(my_read,my_qual), (RC(my_read),my_qual[::-1])]

		nReads += 1
		if nReads%10 == 0:
			print(nReads,'/',nTOT,int(time.time()-tt),'(sec),', int((nTOT-nReads)/(nReads/(time.time()-tt))),'(remaining)')
			#break

		if len(my_read) < MIN_RLEN:
			discard_reads1 += 1
			continue

		all_nice_pairs = []
		for rdi in range(len(rdat_it)):
			readDat      = rdat_it[rdi][0]
			motifs_found = []
			motifs_alts  = []
			if FUZZY:
				# e0
				blacklist = {}
				for i in range(len(MOTIFS)):
					for j in range(0,1+1):
						if MOTIFS[i][j] in readDat:
							newMotifs = [(n.start(0), n.end(0), i, j, 0) for n in re.finditer(MOTIFS[i][j],readDat)]
							motifs_found.extend(newMotifs)
							for n in newMotifs:
								for k in range(n[0],n[1]):
									blacklist[k] = True
				# e1+
				for edit_dist in range(1,FUZZY+1):
					for i in range(len(MOTIFS)):
						for j in range(0,1+1):
							# (start_pos, end_pos, motif_num, GP, edit_dist)
							newMotifs = [(n.start(0), n.end(0), i, j, edit_dist) for n in regex.finditer("("+MOTIFS[i][j]+"){e<=1}", readDat, overlapped=True)]
							#
							for n in newMotifs:
								isInBlack = False
								for k in range(n[0],n[1]):
									if k in blacklist:
										isInBlack = True
										break
								if isInBlack == False:
									motifs_found.append(n)
									for k in range(n[0],n[1]):
										blacklist[k] = True
				#
				motifs_found = sorted(motifs_found + motifs_alts)
				####for n in motifs_found:
				####	print(n)
				####print()
				#
				motifs_found = [(n[0], n[2], n[3]) for n in motifs_found]
			else:
				for i in range(len(MOTIFS)):
					if MOTIFS[i][0] in readDat:
						motifs_found.extend([(n.start(0),i,0) for n in re.finditer(MOTIFS[i][0],readDat)])
					if MOTIFS[i][1] in readDat:
						motifs_found.extend([(n.start(0),i,1) for n in re.finditer(MOTIFS[i][1],readDat)])
				for i in range(len(MOTIF_ALTS)):
					if MOTIF_ALTS[i][2] in readDat:
						motifs_alts.extend([(n.start(0),MOTIF_ALTS[i][0],MOTIF_ALTS[i][1]) for n in re.finditer(MOTIF_ALTS[i][2],readDat)])
			motifs_found = sorted(motifs_found + motifs_alts)

			# prune motifs with low q-scores
			if SKIP_QUALSCORE_FILTERS == False:
				delList = []
				for i in range(len(motifs_found)):
					s1 = motifs_found[i][0]
					e1 = motifs_found[i][0] + len(MOTIFS[motifs_found[i][1]][motifs_found[i][2]])
					if e1 > len(rdat_it[rdi][1]):	# weird edge case where motifs can be found that extend past the actual string end (thanks regex..)
						delList.append(i)
					else:
						q_prio   = [rdat_it[rdi][1][s1+n] for n in MOTIFS[motifs_found[i][1]][4]]
						min_allq = min(rdat_it[rdi][1][s1:e1])
						if len(q_prio):
							min_qpri = min(q_prio)
						else:
							min_qpri = PACBIO_QSCORE_MAX
						if min_allq < QSCORE_MIN_ALL_MOTIF or min_qpri < QSCORE_MIN_PRIO_POS:
							delList.append(i)
						count_min_all[min_allq] += 1
						count_min_pri[min_qpri] += 1
				#print('delList:',delList)
				for i in sorted(delList,reverse=True):
					del motifs_found[i]

			motif = [''.join([symbols[n[2]][n[1]] for n in motifs_found]), [n[0] for n in motifs_found]]

			nice_pairs = []
			for mp in linked_motif_pairs:
				if mp[0] in motif[0] and mp[1] in motif[0]:
					l0 = [i for i in range(len(motif[0])) if motif[0][i] == mp[0]]
					l1 = [i for i in range(len(motif[0])) if motif[0][i] == mp[1]]
					for i in range(len(l0)):
						for j in range(len(l1)):
							ourDist = motif[1][l1[j]] - motif[1][l0[i]]
							cond1   = (mp[2] > 0 and l1[j] > l0[i]) or (mp[2] < 0 and l1[j] < l0[i])
							cond2   = not(any([mmm in motif[0][l0[i]:l1[j]+1] for mmm in mp[3]]))
							cond3   = abs(ourDist-mp[2]) <= DIST_WIGGLE
							#print(mp, (l0[i],l1[j]), (cond1, cond2, cond3), (motif[1][l1[j]], motif[1][l0[i]], ourDist, mp[2], DIST_WIGGLE))
							if cond1 and cond2 and cond3:
								nice_pairs.append([decode[mp[0]], decode[mp[1]], l0[i], l1[j], mp[3]*(mp[3]!=DEAD_CHAR), mp[4]])
			all_nice_pairs.append([(motif[1][n[2]],motif[1][n[3]],n[0],n[1]) for n in nice_pairs])

		# make graphs such that a majority of alignments are FWD
		if len(all_nice_pairs[0]) < len(all_nice_pairs[1]):
			all_nice_pairs = all_nice_pairs[::-1]

		# combine FWD + REV, sort by position
		L = len(my_read)
		sorted_records = sorted([n for n in all_nice_pairs[0]] + [(L-n[1],L-n[0],n[3],n[2]) for n in all_nice_pairs[1]])
		# palindrome detection
		max_pal_so_far = 0
		motif_2_chr    = {}
		if len(all_nice_pairs[0]) and len(all_nice_pairs[1]):
			for n in sorted_records:
				if (n[2],n[3]) not in motif_2_chr and (n[3],n[2]) not in motif_2_chr:
					motif_2_chr[(n[2],n[3])] = chr(len(motif_2_chr)+33)
					motif_2_chr[(n[3],n[2])] = motif_2_chr[(n[2],n[3])]
			motif_str = [motif_2_chr[(n[2],n[3])] for n in sorted_records]
			myPalindromeThresh = int(len(motif_str)*PALINDROME_PERCENTAGE)
			if len(motif_str) >= 4:
				for i in range(myPalindromeThresh,len(motif_str)+1-myPalindromeThresh):
					max_pal_so_far = max([lcs(motif_str[:i][::-1],motif_str[i:]),max_pal_so_far])
					if max_pal_so_far >= myPalindromeThresh:
						break
			#print(my_rname, L, max_pal_so_far, float(max_pal_so_far)/len(motif_str))
			if max_pal_so_far >= myPalindromeThresh:
				discard_reads3 += 1
				continue

		# no questions asked: merge successive nice pairs
		while True:
			foundNone = True
			for j in range(len(sorted_records)):
				for k in range(j+1,len(sorted_records)):
					if sorted_records[j][3] == sorted_records[k][2] and abs(sorted_records[k][0] - sorted_records[j][1]) <= DIST_WIGGLE:
						foundNone = False
						sorted_records[j]  = (sorted_records[j][0], sorted_records[k][1], sorted_records[j][2], sorted_records[k][3])
						del sorted_records[k]
						break
				if not foundNone:
					break
			if foundNone:
				break
		# output
		if len(sorted_records):
			sorted_records = sorted(sorted_records)
			all_data.append([copy.deepcopy(sorted_records),L,my_rname])
		else:
			discard_reads2 += 1

	f.close()

	if len(OUF_READS):
		pickle.dump([PARAMETERS, all_data], open(OUF_READS, "wb"))
	print('discarded reads (<'+str(MIN_RLEN)+'):            ', discard_reads1)
	print('discarded reads (>='+str(MIN_RLEN)+', no motifs):', discard_reads2)
	print('discarded reads (palindromes):     ', discard_reads3)

else:
	print('\nError: Unknown input file type.\n')
	exit(1)


"""//////////////////////////////////////////////////
////////////    CLEANING & FILTERING     ////////////
//////////////////////////////////////////////////"""


print('cleaning reads...')
interesting_connections = {}
cyclic_skips = 0
lowcov_skips = 0
nPlotted = 0
out_clean = []
percentage_explained_hist = []
for i in range(len(all_data)):

	[pairs,L,my_rname] = [copy.deepcopy(all_data[i][0]),all_data[i][1],all_data[i][2]]
	#print(L, my_rname)
	#print('pairs:',pairs)

	P = len(pairs)

	if P == 1:
		[bestScore,bestPath] = [None,[0]]
	else:
		# create graph
		for graph_overlap_this_time in [ALLOWED_OVERLAP_MAX,0]:
			graph_prize    = {}
			graph_weights  = [[None for m in range(P)] for n in range(P)]
			graph_unsorted = []
			gw_dict        = {j:{} for j in range(P)}
			for j in range(P):
				span = abs(int(pairs[j][3].split('_')[1]) - int(pairs[j][2].split('_')[1])) + 1
				graph_prize[j] = float(span)
			for j in range(P):
				tempList = []
				for k in range(P):
					d = pairs[j][0] - pairs[k][1]
					(gp1, gp2) = (pairs[j][2].split('_')[0], pairs[k][2].split('_')[0])
					(s1, e1)   = (int(pairs[j][2].split('_')[1]), int(pairs[j][3].split('_')[1]))
					(s2, e2)   = (int(pairs[k][2].split('_')[1]), int(pairs[k][3].split('_')[1]))
					if j != k and d <= ALLOWED_GAP_MAX and d >= -graph_overlap_this_time:
						if d <= GRAPH_MAX_FREE_GAP and d >= 0:
							graph_weights[j][k] = -0.5
						else:
							if gp1 == gp2 and (e1 < s1) == (e2 < s2):
								graph_weights[j][k] = -1.5
							else:
								graph_weights[j][k] = -2.5
						if d < 0:
							graph_weights[j][k] -= 1.0
						graph_weights[j][k] += graph_prize[k]
						gw_dict[j][k] = graph_weights[j][k]
						tempList.append(k)
				graph_unsorted.append((j,[n for n in tempList]))
			
			# topological sort
			graph_sorted   = []
			graph_unsorted = dict(graph_unsorted)
			graph_fail     = False
			while graph_unsorted:
				acyclic = False
				for node in list(graph_unsorted.keys()):
					edges = graph_unsorted[node]
					for edge in edges:
						if edge in graph_unsorted:
							break
					else:
						acyclic = True
						del graph_unsorted[node]
						graph_sorted.append((node, edges))
				if not acyclic:
					graph_fail = True
					break
			if graph_fail:
				cyclic_skips += 1
			else:
				break
		
		[bestScore,bestPath] = exhaustive_DAG(graph_sorted,graph_weights,graph_prize)
		bestPath = bestPath[::-1]

	# traverse path
	pairs_clean = [pairs[n] for n in bestPath]

	# compute percentage of read explained by observed motifs
	percentage_explained_hist.append(sum([n[1]-n[0]+1 for n in pairs_clean])/float(L))
	#print('percentage_explained:',percentage_explained_hist[-1])
	if percentage_explained_hist[-1] < MIN_READ_COVERAGE:
		lowcov_skips += 1
		#continue

	# clean output strings, should be nice to work with
	out_clean.append([copy.deepcopy(pairs_clean),copy.deepcopy(L),copy.deepcopy(my_rname)])

	# enumerate junctions
	for j in range(1,len(pairs_clean)):
		isFwd1 = (int(pairs_clean[j-1][2].split('_')[1]) < int(pairs_clean[j-1][3].split('_')[1]))
		isFwd2 = (int(pairs_clean[j][2].split('_')[1]) < int(pairs_clean[j][3].split('_')[1]))
		myKey  = (pairs_clean[j-1][3], pairs_clean[j][2], isFwd1, isFwd2)
		myDist = abs(pairs_clean[j][0] - pairs_clean[j-1][1])
		if myKey not in interesting_connections:
			interesting_connections[myKey] = []
		interesting_connections[myKey].append(myDist)

	# dumb strategy
	starts_in = pairs[bestPath[0]][2].split('_')[0]
	ends_in   = pairs[bestPath[-1]][3].split('_')[0]
	starts_in = pairs_clean[0][2].split('_')[0]
	ends_in   = pairs_clean[-1][3].split('_')[0]

	# plot stuff, if desired
	if nPlotted < N_READS_TO_PLOT:
		print('plotting read_'+str(nPlotted+1)+'.png ...')
		my_plotting.plot_read(pairs,pairs_clean,L,my_rname,OUT_DIR+'read_'+str(nPlotted+1)+'.png')
		nPlotted += 1

print('read skipped (unexplained):        ', lowcov_skips)



"""//////////////////////////////////////////////////
////////////     JUNCTION CLUSTERING     ////////////
//////////////////////////////////////////////////"""


print('clustering reads...')

#
sorted_ic = sorted([(len(interesting_connections[k2]),k2) for k2 in interesting_connections.keys()],reverse=True)
ic_clust  = []

for i in range(len(sorted_ic)):
	(me_gp1, me_gp2) = (sorted_ic[i][1][0].split('_')[0], sorted_ic[i][1][1].split('_')[0])
	(me_s, me_e)     = (int(sorted_ic[i][1][0].split('_')[1]), int(sorted_ic[i][1][1].split('_')[1]))
	(me_d1, me_d2)   = (sorted_ic[i][1][2], sorted_ic[i][1][3])
	foundAHome = False
	for j in range(len(ic_clust)):
		a_true_failure  = False
		current_icj_tot = sum([sorted_ic[n][0] for n in ic_clust[j]])
		my_dist_list    = []
		for k in ic_clust[j]:
			(you_gp1, you_gp2) = (sorted_ic[k][1][0].split('_')[0], sorted_ic[k][1][1].split('_')[0])
			(you_s, you_e)     = (int(sorted_ic[k][1][0].split('_')[1]), int(sorted_ic[k][1][1].split('_')[1]))
			(you_d1, you_d2)   = (sorted_ic[k][1][2], sorted_ic[k][1][3])
			if me_gp1 != you_gp1 or me_gp2 != you_gp2 or me_d1 != you_d1 or me_d2 != you_d2:
				a_true_failure = True
				break
			my_dist_list.append((abs(me_s-you_s)+abs(me_e-you_e), sorted_ic[k][0]))
		yay = sum([my_dist_list[n][1] for n in range(len(my_dist_list)) if my_dist_list[n][0] <= MAX_SHIFT_DIST])
		nay = sum([my_dist_list[n][1] for n in range(len(my_dist_list)) if my_dist_list[n][0] > MAX_SHIFT_DIST])
		if a_true_failure == False and yay/float(current_icj_tot) >= MAX_SHIFT_DIST_FRAC:
			ic_clust[j].append(i)
			foundAHome = True

	if not foundAHome:
		ic_clust.append([i])

#
happy_families = {}
span_by_fam    = {}
dist_by_fam    = {}
ic_byCount = sorted([(sum([sorted_ic[n][0] for n in ic_clust[m]]),m) for m in range(len(ic_clust))],reverse=True)
for i in range(len(ic_byCount)):
	spanList = [[],[]]
	for j in ic_clust[ic_byCount[i][1]]:
		spanList[0].append(int(sorted_ic[j][1][0].split('_')[1]))
		spanList[1].append(int(sorted_ic[j][1][1].split('_')[1]))
		(mygp1, mygp2) = (sorted_ic[j][1][0][0], sorted_ic[j][1][1][0])
		(myfr1, myfr2) = (sorted_ic[j][1][2], sorted_ic[j][1][3])
		happy_families[sorted_ic[j][1]] = i
	span_by_fam[i] = [mygp1, myfr1, min(spanList[0]), max(spanList[0]), mygp2, myfr2, min(spanList[1]), max(spanList[1])]
	dist_by_fam[i] = [l_2_cd(spanList[0]), l_2_cd(spanList[1])]

#
all_fam = {}
for i in range(len(out_clean)):
	[pairs_clean,L,my_rname] = out_clean[i]
	[myS, myE]               = [pairs_clean[0][2], pairs_clean[-1][3]]
	all_keys                 = []
	for j in range(1,len(pairs_clean)):
		isFwd1 = (int(pairs_clean[j-1][2].split('_')[1]) < int(pairs_clean[j-1][3].split('_')[1]))
		isFwd2 = (int(pairs_clean[j][2].split('_')[1]) < int(pairs_clean[j][3].split('_')[1]))
		myKey  = (pairs_clean[j-1][3], pairs_clean[j][2], isFwd1, isFwd2)
		myDist = abs(pairs_clean[j][0] - pairs_clean[j-1][1])
		all_keys.append(happy_families[myKey])
	if tuple(all_keys) not in all_fam:
		all_fam[tuple(all_keys)] = []
	all_fam[tuple(all_keys)].append([i,myS,myE])


"""//////////////////////////////////////////////////////////////////
////////////     PLOTTING & STRUCTURE IDENTIFICATION     ////////////
//////////////////////////////////////////////////////////////////"""


af_byCount  = sorted([[len(all_fam[k]),k] for k in all_fam.keys()],reverse=True)
unaccounted = 0
nPlotted    = 0
af_tot      = float(sum([n[0] for n in af_byCount]))
rolling_tot = af_tot
form_names  = {1:'A',2:'B',3:'C',4:'D',5:'A-like',6:'B-like',7:'C-like',8:'D-like',9:'other'}
form_names  = {1:'normal gene',2:'gene/pseudogene chimera',3:'pseudogene/gene chimera',4:'normal pseudogene',5:'normal gene (weak)',6:'gene/pseudogene chimera (weak)',7:'pseudogene/gene chimera (weak)',8:'normal pseudogene (weak)',9:'other'}
form_counts = {n:0 for n in form_names.keys()}
form_counts_weak = {n:0 for n in form_names.keys()}
forms_max_strong = {1:2, 2:2, 3:2, 4:2, 5:4, 6:2, 7:2, 8:4, 9:999}
vcf_dat     = []
rname_dat   = []
REGION_NAME_DICT = {('G',True):'CYP21A2_FWD', ('G',False):'CYP21A2_REV', ('P',True):'CYP21A1P_FWD', ('P',False):'CYP21A1P_REV'}
for i in range(len(af_byCount)):
	if af_byCount[i][0] < MIN_READS_PER_CLUSTER:
		unaccounted += af_byCount[i][0]

	# singletons
	if af_byCount[i][1] == ():
		#print(af_byCount[i])

		allReads      = [out_clean[m[0]] for m in all_fam[af_byCount[i][1]] if m[1][0] == 'G']
		allReads_weak = []
		for j in range(len(allReads)-1,-1,-1):
			inds = int(allReads[j][0][0][2].split('_')[1])
			indf = int(allReads[j][0][0][3].split('_')[1])
			if min([inds,indf]) > HIGH_CONF_ANCHOR or max([inds,indf]) < len(MOTIFS)-HIGH_CONF_ANCHOR-1:
				allReads_weak.append(copy.deepcopy(allReads[j]))
				del allReads[j]

		my_form                    = 1
		form_counts[my_form]      += len(allReads)
		form_counts_weak[my_form] += len(allReads_weak)
		strong_percent             = 100.*(len(allReads)/af_tot)
		weak_percent               = 100.*(len(allReads_weak)/af_tot)

		if len(allReads)+len(allReads_weak)  >= MIN_READS_PER_CLUSTER:
			# plotting
			print('plotting cluster_'+str(nPlotted+1)+'.png ...')
			if EXTRA_PLOTS:
				my_plotting.plot_AD_hist(allReads,allReads_weak,plotTitle='gene motif coverage',saveName=OUT_DIR+'hist_G.png')
			if len(allReads) > 0:
				myTitle = 'reads: '+str(len(allReads))+' ({0:0.2f}%)'.format(strong_percent)+', '+form_names[my_form]
				my_plotting.plot_multiple(allReads,plotTitle=myTitle,saveName=OUT_DIR+'cluster_'+str(nPlotted+1)+'.png')
			if EXTRA_PLOTS and len(allReads_weak) > 0:
				myTitle = 'reads: '+str(len(allReads_weak))+' ({0:0.2f}%)'.format(weak_percent)+', '+form_names[my_form]+' (weak)'
				my_plotting.plot_multiple(allReads_weak,plotTitle=myTitle,saveName=OUT_DIR+'cluster_weak_'+str(nPlotted+1)+'.png')
			# output vcf dat
			vcf_dat.append([nPlotted, form_names[my_form], len(allReads)+len(allReads_weak), strong_percent+weak_percent, [], []])
			# output names of reads that strongly support this cluster for subsequent genotyping
			if len(allReads) > 0:
				rname_dat.append(('cluster_'+str(nPlotted+1), form_names[my_form], [], [], [n[2] for n in allReads]))
			nPlotted += 1

		allReads      = [out_clean[m[0]] for m in all_fam[af_byCount[i][1]] if m[1][0] == 'P']
		allReads_weak = []
		for j in range(len(allReads)-1,-1,-1):
			inds = int(allReads[j][0][0][2].split('_')[1])
			indf = int(allReads[j][0][0][3].split('_')[1])
			if min([inds,indf]) > HIGH_CONF_ANCHOR or max([inds,indf]) < len(MOTIFS)-HIGH_CONF_ANCHOR-1:
				allReads_weak.append(copy.deepcopy(allReads[j]))
				del allReads[j]

		my_form                    = 4
		form_counts[my_form]      += len(allReads)
		form_counts_weak[my_form] += len(allReads_weak)
		strong_percent             = 100.*(len(allReads)/af_tot)
		weak_percent               = 100.*(len(allReads_weak)/af_tot)

		if len(allReads)+len(allReads_weak) >= MIN_READS_PER_CLUSTER:
			print('plotting cluster_'+str(nPlotted+1)+'.png ...')
			# plotting
			if EXTRA_PLOTS:
				my_plotting.plot_AD_hist(allReads,allReads_weak,plotTitle='gene motif coverage',saveName=OUT_DIR+'hist_P.png')
			if len(allReads) > 0:
				myTitle = 'reads: '+str(len(allReads))+' ({0:0.2f}%)'.format(strong_percent)+', '+form_names[my_form]
				my_plotting.plot_multiple(allReads,plotTitle=myTitle,saveName=OUT_DIR+'cluster_'+str(nPlotted+1)+'.png')
			if EXTRA_PLOTS and len(allReads_weak) > 0:
				myTitle = 'reads: '+str(len(allReads_weak))+' ({0:0.2f}%)'.format(weak_percent)+', '+form_names[my_form]+' (weak)'
				my_plotting.plot_multiple(allReads_weak,plotTitle=myTitle,saveName=OUT_DIR+'cluster_weak_'+str(nPlotted+1)+'.png')
			# output vcf dat
			vcf_dat.append([nPlotted, form_names[my_form], len(allReads)+len(allReads_weak), strong_percent+weak_percent, [], []])
			# output names of reads that strongly support this cluster for subsequent genotyping
			if len(allReads) > 0:
				rname_dat.append(('cluster_'+str(nPlotted+1), form_names[my_form], [], [], [n[2] for n in allReads]))
			nPlotted += 1

		rolling_tot -= af_byCount[i][0]
		#print('rolling:', int(rolling_tot), '({0:0.2f}%)'.format(100.*(rolling_tot/af_tot)))
	else:
		myConn = [span_by_fam[n] for n in af_byCount[i][1]]
		bpDist = [dist_by_fam[n] for n in af_byCount[i][1]]

		# parse breakpoint/region info for vcf and rname output
		my_vcf_dat_regions = []
		my_vcf_dat_breaks  = []
		for j in range(len(myConn)):
			if j == 0:
				my_vcf_dat_regions.append(REGION_NAME_DICT[(myConn[j][0],myConn[j][1])])
			my_vcf_dat_regions.append(REGION_NAME_DICT[(myConn[j][4],myConn[j][5])])
			# simple simple simple logic, just use mode! todo: investigate more sophisticated choice
			bp1m = bpDist[j][0][0][0]
			bp2m = bpDist[j][1][0][0]
			my_vcf_dat_breaks.extend([MOTIFS[bp1m][2+1*(myConn[j][0] == 'P')], MOTIFS[bp2m][2+1*(myConn[j][4] == 'P')]])

		# determine label
		my_form = 9
		if len(myConn) == 1:
			if myConn[0][1] and myConn[0][5]:					# is fwd
				if min(myConn[0][2:4]) < max(myConn[0][6:8]):	# makes progress along gene
					if myConn[0][0] == 'G' and myConn[0][4] == 'G': my_form = 5	# imperfect A
					if myConn[0][0] == 'G' and myConn[0][4] == 'P': my_form = 2	# B
					if myConn[0][0] == 'P' and myConn[0][4] == 'G': my_form = 3	# C
					if myConn[0][0] == 'P' and myConn[0][4] == 'P': my_form = 8	# imperfect D
		else:
			if len(myConn) >= 4:
				my_form = 9	# too much noise, call it "other"
			else:
				if myConn[0][1] and myConn[-1][5]:
					if min(myConn[0][2:4]) < max(myConn[-1][6:8]):
						if myConn[0][0] == 'G' and myConn[-1][4] == 'G': my_form = 5	# imperfect A
						if myConn[0][0] == 'G' and myConn[-1][4] == 'P': my_form = 6	# imperfect B
						if myConn[0][0] == 'P' and myConn[-1][4] == 'G': my_form = 7	# imperfect C
						if myConn[0][0] == 'P' and myConn[-1][4] == 'P': my_form = 8	# imperfect D

		# penalize >2 different B/C forms when performing high-conf summary counts
		tooManyDifferent = (forms_max_strong[my_form] <= 0)
		if my_form in [1,5]:
			forms_max_strong[1] -= 1
			forms_max_strong[5] -= 1
		elif my_form in [2,6]:
			forms_max_strong[2] -= 1
			forms_max_strong[6] -= 1
		elif my_form in [3,7]:
			forms_max_strong[3] -= 1
			forms_max_strong[7] -= 1
		elif my_form in [4,8]:
			forms_max_strong[4] -= 1
			forms_max_strong[8] -= 1

		if tooManyDifferent:	# if we've seen too many of this form already, call them low confidence
			allReads      = []
			allReads_weak = [out_clean[m[0]] for m in all_fam[af_byCount[i][1]]]
		else:
			allReads      = [out_clean[m[0]] for m in all_fam[af_byCount[i][1]]]
			allReads_weak = []
			for j in range(len(allReads)-1,-1,-1):
				inds1 = int(allReads[j][0][0][2].split('_')[1])
				indf1 = int(allReads[j][0][0][3].split('_')[1])
				inds2 = int(allReads[j][0][-1][2].split('_')[1])
				indf2 = int(allReads[j][0][-1][3].split('_')[1])
				if min([inds1,indf1]) > HIGH_CONF_ANCHOR or max([inds2,indf2]) < len(MOTIFS)-HIGH_CONF_ANCHOR-1:
					allReads_weak.append(copy.deepcopy(allReads[j]))
					del allReads[j]

		form_counts[my_form]      += len(allReads)
		form_counts_weak[my_form] += len(allReads_weak)
		strong_percent = 100.*(len(allReads)/af_tot)
		weak_percent   = 100.*(len(allReads_weak)/af_tot)

		if len(allReads)+len(allReads_weak) >= MIN_READS_PER_CLUSTER:
			print('plotting cluster_'+str(nPlotted+1)+'.png ...')
			# plotting
			if len(allReads) > 0:
				myTitle = 'reads: '+str(len(allReads))+' ({0:0.2f}%)'.format(strong_percent)+', '+form_names[my_form]
				my_plotting.plot_multiple(allReads,plotTitle=myTitle,saveName=OUT_DIR+'cluster_'+str(nPlotted+1)+'.png')
			if EXTRA_PLOTS and len(allReads_weak) > 0:
				myTitle = 'reads: '+str(len(allReads_weak))+' ({0:0.2f}%)'.format(weak_percent)+', '+form_names[my_form]+' (weak)'
				my_plotting.plot_multiple(allReads_weak,plotTitle=myTitle,saveName=OUT_DIR+'cluster_weak_'+str(nPlotted+1)+'.png')
			# output vcf dat
			vcf_dat.append([nPlotted, form_names[my_form], len(allReads)+len(allReads_weak), strong_percent+weak_percent, my_vcf_dat_breaks, my_vcf_dat_regions])
			# output names of reads that strongly support this cluster for subsequent genotyping
			if len(allReads) > 0:
				rname_dat.append(('cluster_'+str(nPlotted+1), form_names[my_form], my_vcf_dat_breaks, my_vcf_dat_regions, [n[2] for n in allReads]))
			nPlotted += 1

		rolling_tot -= af_byCount[i][0]
		#print('rolling:', int(rolling_tot), '({0:0.2f}%)'.format(100.*(rolling_tot/af_tot)))
print('reads not plotted:', unaccounted)

print('writing output VCF...')
writeVCF(vcf_dat, OUT_DIR+'summary.vcf')

print('writing output readnames (sorted by cluster)...')
f = open(OUT_DIR+'readClusters.txt','w')
for i in range(len(rname_dat)):
	f.write(rname_dat[i][0] + '\t' + rname_dat[i][1] + '\t')
	f.write(' '.join([str(n) for n in rname_dat[i][2]]) + '\t')
	f.write(' '.join([str(n) for n in rname_dat[i][3]]) + '\t')
	f.write(' '.join([str(n) for n in rname_dat[i][4]]) + '\n')
f.close()


"""//////////////////////////////////////////////////
////////////        PRINT RESULTS        ////////////
//////////////////////////////////////////////////"""

BAR = '============================'
SPC = '                            \n'

print(SPC + 'RESULTS (HIGH-CONF):        \n' + SPC)
if sum(list(form_counts.values())) > 0:
	print_nicely({form_names[k]:form_counts[k] for k in form_counts.keys()},
		         order=[form_names[k] for k in sorted(form_counts.keys())],
		         end_boof=len(BAR),
		         custom_tot=sum(list(form_counts.values())+list(form_counts_weak.values())))

print(SPC + 'RESULTS (LOW-CONF):         \n' + SPC)
if sum(list(form_counts_weak.values())) > 0:
	print_nicely({form_names[k]:form_counts_weak[k] for k in form_counts_weak.keys()},
		         order=[form_names[k] for k in sorted(form_counts_weak.keys())],
		         end_boof=len(BAR),
		         custom_tot=sum(list(form_counts.values())+list(form_counts_weak.values())))

summarized = {'A form':form_counts[1]+form_counts[5],
              'B form':form_counts[2]+form_counts[6],
              'C form':form_counts[3]+form_counts[7],
              'D form':form_counts[4]+form_counts[8]}

print(SPC + 'SUMMARY (HIGH-CONF):        \n' + SPC)
if sum(list(summarized.values())) > 0:
	print_nicely(summarized,end_boof=len(BAR))
othstr = 'other:   '+str(form_counts[9])#,' [{0:0.2f}]'.format(form_counts[9]/float(sum([form_counts[n] for n in range(1,9+1)])))
print(othstr+' '*(max([0,len(BAR)-len(othstr)])))

summarized = {'A form':form_counts[1]+form_counts[5]+form_counts_weak[1]+form_counts_weak[5],
              'B form':form_counts[2]+form_counts[6]+form_counts_weak[2]+form_counts_weak[6],
              'C form':form_counts[3]+form_counts[7]+form_counts_weak[3]+form_counts_weak[7],
              'D form':form_counts[4]+form_counts[8]+form_counts_weak[4]+form_counts_weak[8]}

print(SPC + 'SUMMARY (ALL):              \n' + SPC)
if sum(list(summarized.values())) > 0:
	print_nicely(summarized,end_boof=len(BAR))
othstr = 'other:   '+str(form_counts[9]+form_counts_weak[9])
print(othstr+' '*(max([0,len(BAR)-len(othstr)])))
print(SPC)

