#
#	construct MOTIFS data structure, using manually derived motif sequences and motif positions,
#	as well as sequences related to flanking primers and what-not.
#
#	all this code is very specific to PMS2. It would need to be
#	extensively re-derived if you wish to apply this methodology to another gene/pseudogene
#
#	script usage: python pbcyp_makeMotifs.py output.p
#

import sys
import cPickle as pickle

# main data structure containing motif data
#
# (gene motif, pseudogene motif, gene coordinates, pseudogene coordinates, high-priority motif coordinates, gene alts, pseudogene alts)
#
# coordinates don't correspond to reference position, due to strand directions being annoying
MOTIFS = [('AAAACAAACAGTAA',   'AAAAAAGCACTCAGTAA', 29428, 9902,  []),
          ('GAGAGGGATGG',      'GAGAAGGATGG',       29453, 9930,  []),
          ('CAGGCAGAGGA',      'CAGGGAGAGGA',       29474, 9951,  []),
          ('TAAAAGCTGTC',      'TAAATGCTGTC',       29503, 9980,  []),
          ('TTCTGTATTAG',      'TTCTCTATTAG',       29564, 10041, []),
          ('ACTTGAACTTTTTT',   'ACTTTAACTTTTTTT',   29687, 10164, []),
          ('TTTTATGGGTG',      'TTTTGTGGGTG',       29866, 10344, []),
          ('AGAGTGTTTCG',      'AGAGGGTTTCG',       30075, 10558, []),
          ('TCGCTGTGTTG',      'TCGCCGTGTTG',       30083, 10566, []),
          ('AGGCCGGTCTCAA',    'AGGCTGATCTCAA',     30097, 10580, []),
          ('CCTCGGGCTCCTA',    'CCTCAGCCTCCTA',     30136, 10619, []),
          ('GCTACTTTAAC',      'GCTATTTTAAC',       30192, 10675, []),
          ('CTGAGTAGTGTTGT',   'CTGAATTATGTTGT',    30269, 10752, []),
          ('TATTTGTGTTG',      'TATTCGTGTTG',       30283, 10766, []),
          ('GCTTTGATTTT',      'GCTTCGATTTT',       30478, 10961, []),
          ('GCCTGGCTAAC',      'GCCTCGCTAAC',       30617, 11100, []),
          ('AAATCAGAACT',      'AAATTAGAACT',       30825, 11310, []),
          ('AGACCAGTTTTTT',    'AGACTGTTTTTTT',     30883, 11368, []),
          ('AAAAAGTTCCTCCAAA', 'AAAAAGGTTCCCCCAAA', 30914, 11399, []),
          ('TGTGAGTTTCCTC',    'TGTGACTGTCCTC',     30938, 11424, []),
          ('TTCAGGTGTTC',      'TTCAGATGTTC',       30954, 11440, []),
          ('GGAGCGCAGTG',      'GGAGCACAGTG',       31010, 11503, []),
          ('GGTCTGTATCT',      'GGTCTCTATCT',       31179, 11671, []),
          ('TAAAAACTGGA',      'TAAAAGCTGGA',       31386, 11878, []),
          ('GTCTGCGTTCT',      'GTCTGTGTTCT',       31677, 12169, []),
          ('CGGAACTTGCA',      'CGGAAGTTGCA',       31959, 12451, []),
          ('CTGATGTCCAA',      'CTGATTTCCAA',       32590, 13081, []),
          ('CGCTGAGGAGA',      'CGCTGGGGAGA',       33833, 14322, []),
          ('CGCCTGGTGC',       'CGCCTGGGTGC',       34696, 15185, []),
          ('AAATGAACCT',       'AAATGAAACCT',       35787, 16277, [])]

# print out motifs in bed format
#for i in range(len(MOTIFS)):
#	print 'chr7\t'+str(MOTIFS[i][3])+'\t'+str(MOTIFS[i][3]+len(MOTIFS[i][1]))+'\tmotif_'+str(i)+'_'+MOTIFS[i][1]+'_('+MOTIFS[i][0]+')'
#for i in range(len(MOTIFS)):
#	print 'chr7\t'+str(MOTIFS[i][2])+'\t'+str(MOTIFS[i][2]+len(MOTIFS[i][0]))+'\tmotif_'+str(i)+'_'+MOTIFS[i][0]+'_('+MOTIFS[i][1]+')'
#exit(1)

# motifs that can be altered via common snps
MOTIF_ALTS = []

N_PRIMER_MOTIFS = 0

# trim motifs for more hits but less specificity, if desired
MOTIF_TRIM = 1 # flank = 5 - MOTIF_TRIM
if MOTIF_TRIM > 0:
	for i in range(len(MOTIFS)):
		for k in range(len(MOTIFS[i][4])-1,-1,-1):
			if MOTIFS[i][4][k] >= min([len(MOTIFS[i][0]),len(MOTIFS[i][1])]) - MOTIF_TRIM:
				del MOTIFS[i][4][k]
		MOTIFS[i] = (MOTIFS[i][0][MOTIF_TRIM:-MOTIF_TRIM], MOTIFS[i][1][MOTIF_TRIM:-MOTIF_TRIM], MOTIFS[i][2]+MOTIF_TRIM, MOTIFS[i][3]+MOTIF_TRIM, [n-MOTIF_TRIM for n in MOTIFS[i][4] if n-MOTIF_TRIM >= 0])
	for i in range(len(MOTIF_ALTS)):
		MOTIF_ALTS[i] = (MOTIF_ALTS[i][0], MOTIF_ALTS[i][1], MOTIF_ALTS[i][2][MOTIF_TRIM:-MOTIF_TRIM])

# PMS2-specific metadata for plotting
# (hg38 coords)
MY_EXONS = [(5970924, 5973542), (5977587, 5977757), (5978595, 5978696), (5982823, 5982991),
             (5986758, 5987620), (5989799, 5989955), (5991972, 5992057), (5995533, 5995633),
             (5997325, 5997423), (5999107, 5999275), (6002452, 6002636), (6003689, 6003792),
             (6003971, 6004058), (6005891, 6006031), (6008996, 6009049)]
MY_EXONS = []
MY_XLAB  = {'G':'PMS2', 'P':'PMS2CL'}
PLOTTING_META = {'exons':MY_EXONS, 'xlabel':MY_XLAB}

#
# save output
#
pickle.dump([MOTIFS, MOTIF_ALTS, N_PRIMER_MOTIFS, MOTIF_TRIM, PLOTTING_META], open(sys.argv[1], "wb"))


