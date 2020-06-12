import datetime

VCF_HEADER = ['##fileformat=VCFv4.2',
			  '##source=pacbio_cyp21a2',
			  '##fileDate=99999999',
			  '##contig=<ID=chr6,length=170805979>',
			  '##ALT=<ID=A,Description="A form">',
			  '##ALT=<ID=B,Description="B form">',
			  '##ALT=<ID=C,Description="C form">',
			  '##ALT=<ID=D,Description="D form">',
			  '##ALT=<ID=A-like,Description="starts and ends in CYP21A2">',
			  '##ALT=<ID=B-like,Description="starts in CYP21A2, ends in CYP21A1P">',
			  '##ALT=<ID=C-like,Description="starts in CYP21A1P, ends in CYP21A2">',
			  '##ALT=<ID=D-like,Description="starts and ends in CYP21A1P">',
			  '##INFO=<ID=BREAKPOINTS,Number=.,Type=Integer,Description="Highest confidence breakpoint positions">',
			  '##INFO=<ID=REGIONS,Number=.,Type=String,Description="Indicates what portion of event is in gene or pseudogene">',
			  '##INFO=<ID=ZMW,Number=1,Type=Integer,Description="Number of ZMWs supporting SV.">',
			  '##INFO=<ID=ZMW_PERCENT,Number=1,Type=Integer,Description="Percentage of total ZMWs supporting SV.">',
			  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmySample']

# starting coordinate for all rearrangements arbitrarily chosen
# to be the position of the U-form primer in CYP21A2 region:
COORD_CHR   = 'chr6'
COORD_START = 32038152

def writeVCF(records,out_fileName):
	"""
	write output VCF file

	records - [cluster_num, form_label, nReads, percentReads, breakpointList, regionList]
	"""

	VCF_HEADER[2] = '##fileDate=' + datetime.datetime.now().strftime("%y%m%d")
	f = open(out_fileName,'w')
	for line in VCF_HEADER:
		f.write(line + '\n')
	for n in records:
		dat = [COORD_CHR, str(COORD_START), str(n[0]), 'N', '<'+n[1]+'>', '.', 'PASS']
		myInfo = 'ZWM='+str(n[2])+';'+'ZMW_PERCENT={0:0.2f}%;'.format(n[3])
		if len(n[4]):
			myInfo += 'BREAKPOINTS=' + ','.join([str(m) for m in n[4]]) + ';'
		if len(n[5]):
			myInfo += 'REGIONS=' + ','.join(n[5]) + ';'
		dat.append(myInfo)
		f.write('\t'.join(dat)+'\t'+'GT'+'\t'+'./.'+'\n')
	f.close()
