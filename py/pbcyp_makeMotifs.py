#
#	construct MOTIFS data structure, using manually derived motif sequences and motif positions,
#	as well as sequences related to flanking primers and what-not.
#
#	all this code is very specific to CYP21A2, and our probe design. It would need to be
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
# coordinates are hg38::chr6
MOTIFS = [('CTTCTTGATGG',           'CTTCTCGATGG',           32038207, 32005472, [5]                                   ),
		  ('GTGATCAATTTTTTTGAAAT',  'GTGATTAATTTTTTTTGAAAT', 32038218, 32005483, [5,15],  ['GTGATCAACTTTTTTGAAAT'], [] ),
		  ('GTGGGCGGGTC',           'GTGGGTGGGTC',           32038291, 32005557, [5]                                   ),
		  ('TGGGAGGGTACCTG',        'TGGGAAGGCACCTG',        32038304, 32005570, [5,8]                                 ),
		  ('CCTGAAGGTGG',           'CCTGAGGGTGG',           32038314, 32005580, [5]                                   ),
		  ('CGTCTCGCCAT',           'CGTCTTGCCAT',           32038413, 32005679, [5]                                   ),
		  ('CCTCCCGCCTC',           'CCTCCTGCCTC',           32038508, 32005774, [5]                                   ),
		  ('CTTAAAAAAATTT',         'CTTAAACAAATTT',         32038837, 32006103, [6],     [], ['CTTAAACAATTTT']        ),
		  ('TTTTTAAGAGA',           'TTTTTGTTAGAGA',         32038851, 32006118, [5,6,7]                               ),
		  ('ATGGGTTCTTG',           'ATGGGGTCTTG',           32038861, 32006130, [5]                                   ),
		  ('CTATGCTGCCC',           'CTATGTTGCCC',           32038872, 32006141, [5]                                   ),
		  ('GTCTTAAATTC',           'GTCTTGAATTC',           32038889, 32006158, [5]                                   ),
		  ('TTCCTAGTCTC',           'TTCCTGGTCTC',           32038897, 32006166, [5]                                   ),
		  ('CTCAAATGATC',           'CTCAAGTGATC',           32038905, 32006174, [5]                                   ),
		  ('ACCTCAGCCTC',           'ACCTCGGCCTC',           32038921, 32006190, [5]                                   ),
		  ('AAGTGTGAGCC',           'AAGTGGGAGCC',           32038932, 32006201, [5]                                   ),
		  ('ACCTTTGGGGC',           'ACCTTCGGGGG',           32038943, 32006212, [5]                                   ),
		  ('GGGGCATCCCC',           'GGGGCTTCCCC',           32038949, 32006219, [5],     [], ['GGGGCTTCTCC']          ),
		  #('AATCCAGGTC',            'AATCCTCCAGGTC',         32038960, 32006230, [5,6,7], [], ['AATCCTTCAGGTC']        ),	# [32007021, 32038960, 32039757] [32006230]
		  ('AGGTCCCTGGA',           'AGGTCACTGGA',           32038965, 32006238, [5]                                   ),
		  ('CTCTTGGGGGGGCATAT',     'CTCTTGGGGGGCATAT',      32038978, 32006251, [11]                                  ),
		  ('TATCTGGTGGGGAGAAAGCA',  'TATCTTCAGGAGAAGAAGCA',  32038992, 32006264, [5,6,7,8]                             ),
		  ('GCAGGGGTTGGGGAGG',      'GCAGGTGTTGAGGAGG',      32039009, 32006281, [5,10]                                ),
		  ('GAGGCCGAAGA',           'GAGGCAGAAGA',           32039021, 32006293, [5]                                   ),
		  ('CCCTCAGCTGCCTTC',       'CCCTCGGCTTCCTTG',       32039040, 32006312, [5,9]                                 ),
		  ('GCCTTCATCAGT',          'TCCTTGGTCAGT',          32039049, 32006321, [0,5,6]                               ),
		  ('CCCCACCTCCT',           'CCCCAGCTCCT',           32039075, 32006347, [5],     ['CCCCACTTCCT'], []          ),
		  ('GTCTAGGAACT',           'GTCTAAGAACT',           32039103, 32006375, [5]                                   ),
		  ('CTGTCCTTGGG',           'CTGTCGTTGGT',           32039122, 32006394, [5],     ['CTTTCCTTGGG'], []          ),
		  ('CTTGGGAGACTACTCCCTGCT', 'GTTGGTCTCTGCT',         32039127, 32006399, [5,6,7,8]                             ),
		  ('CATCATCTGTT',           'CATCAACTGTT',           32039420, 32006684, [5]                                   ),
		  ('ACTCCCTCCTT',           'ACTCCATCCTT',           32039525, 32006789, [5]                                   ),
		  ('CCTTTTCTGGC',           'CCTTTCCTGGC',           32039532, 32006796, [5]                                   ),
		  ('CAGGACGACAA',           'CAGGAGGACAA',           32039542, 32006806, [5]                                   ),
		  ('AGGGATCACATCGTGGAGA',   'AGGGACCACAACGAGGAGA',   32039796, 32007060, [5,10,13]                             ),
		  ('GGAGATGCAGC',           'GGAGAAGCAGC',           32039810, 32007074, [5]                                   ),
		  ('ACTGTACGTGGA',          'ACTGTGTGTGGA',          32039841, 32007105, [5,6]                                 ),
		  ('AGCCTCGTGGC',           'AGCCTGGTGGC',           32040007, 32007271, [5]                                   ),
		  ('GGCACGTGCAC',           'GGCACTTGCAC',           32040104, 32007368, [5]                                   ),
		  ('CGTGGTTTTTTTGCTTC',     'CGTGGTTTTTTTTGCTTC',    32040177, 32007441, [4,5,11,12]                           ),
		  ('TCCTGGGGACA',           'TCCTGCGGACA',           32040210, 32007475, [5]                                   ),
		  ('GACTGCAGGAG',           'GACTGTAGGAG',           32040415, 32007680, [5]                                   ),
		  ('GCCTGCGGCCC',           'GCCTGTGGCCC',           32040529, 32007794, [5]                                   ),
		  ('CTGCCGTGAAA',           'CTGCCATGAAA',           32040791, 32008056, [5]                                   ),
		  ('AGCCCGGGCCA',           'AGCCCAGGCCA',           32041113, 32008378, [5]                                   ),
		  ('CCAGAGCCAGT',           'CCAGAACCAGT',           32041121, 32008386, [5]                                   ),
		  ('TGCTCTTCCCG',           'TGCTCCTCCCG',           32041496, 32008761, [5],     [], ['TGCTCCTCCCA']          ),
		  ('GAGGTAGCTCC',           'GAGGTGGCTCC',           32041518, 32008783, [5]                                   ),
		  ('CCCTCCGCTGCAGA',        'CCCTCTGCCGCAGA',        32041568, 32008833, [5,8]                                 ),
		  ('TTAATTCTGAG',           'TTAATCCTGAG',           32041592, 32008857, [5]                                   ),
		  ('GCTGGCCCTTT',           'GCTGGTCCTTT',           32041602, 32008867, [5]                                   )]
		  #('CAGCCACTCAG',           'CAGCCGCTCAG',           32041673, 32008938, [5]                                   )]	# [32041673] [32007234, 32008938, 32039970]

# print out motifs in bed format
#for i in range(len(MOTIFS)):
#	print 'chr6\t'+str(MOTIFS[i][3])+'\t'+str(MOTIFS[i][3]+len(MOTIFS[i][1]))+'\tmotif_'+str(i)+'_'+MOTIFS[i][1]+'_('+MOTIFS[i][0]+')'
#for i in range(len(MOTIFS)):
#	print 'chr6\t'+str(MOTIFS[i][2])+'\t'+str(MOTIFS[i][2]+len(MOTIFS[i][0]))+'\tmotif_'+str(i)+'_'+MOTIFS[i][0]+'_('+MOTIFS[i][1]+')'
#exit(1)

# motifs that can be altered via common snps
MOTIF_ALTS = []
for i in range(len(MOTIFS)):
	if len(MOTIFS[i]) >= 7:
		for j in range(len(MOTIFS[i][5])):
			MOTIF_ALTS.append((i,0,MOTIFS[i][5][j]))
		for j in range(len(MOTIFS[i][6])):
			MOTIF_ALTS.append((i,1,MOTIFS[i][6][j]))

# append some motifs to the end representing primers
# (this is to capture novel motifs at the start/end of gene/pseudogene)
PRIMER_FWD = [('CAGAAAGCTGAC', 'CAGAAAGCTGAC', 32038152, 32005417, [5]),
			  ('AAGCTGACTCTG', 'AAGCTGACTCTG', 32038156, 32005421, [5]),
			  ('TGACTCTGGATG', 'TGACTCTGGATG', 32038160, 32005425, [5]),
			  ('TCTGGATGCAGG', 'TCTGGATGCAGG', 32038164, 32005429, [5])]
PRIMER_REV = [('GTTGAGGTTGGC', 'GTTGAGGTTGGC', 32041800, 32009065, [5]),
			  ('AGGTTGGCGTAG', 'AGGTTGGCGTAG', 32041804, 32009069, [5]),
			  ('TGGCGTAGTGGC', 'TGGCGTAGTGGC', 32041808, 32009073, [5]),
			  ('GTAGTGGCAGTT', 'GTAGTGGCAGTT', 32041812, 32009077, [5])]

N_PRIMER_MOTIFS = len(PRIMER_FWD)

MOTIFS = PRIMER_FWD + MOTIFS + PRIMER_REV

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

# CYP21A2-specific metadata for plotting
# (hg38 coords)
MY_EXONS = [(32038315,32038624), (32038721,32038811), (32039093,32039248),
            (32039355,32039457), (32039545,32039647), (32039748,32039835),
            (32040004,32040205), (32040405,32040584), (32040667,32040771),
            (32040868,32041670)]
MY_XLAB  = 'CYP21A1P <----            ----> CYP21A2 '
PLOTTING_META = {'exons':MY_EXONS, 'xlabel':MY_XLAB}

#
# save output
#
pickle.dump([MOTIFS, MOTIF_ALTS, N_PRIMER_MOTIFS, MOTIF_TRIM, PLOTTING_META], open(sys.argv[1], "wb"))


