import numpy as np

import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def getColor(i,N,colormap='jet'):
	"""
	grab a color from a colormap
	
	i - index of array in colormap array
	N - total length of colormap array
	returns: 4-tuple: (r,g,b,a)
	"""
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i)
	return colorVal


class plotting_object:
	"""
	class for handling all the plot generation
	"""
	def __init__(self, MOTIFS, PLOTTING_META):
		self.L_MOTIF  = len(MOTIFS)
		#self.BOOF_GP  = 1000
		#self.BOOF_XE  = 100
		#self.X_COORD  = [n[3]-MOTIFS[0][3] for n in MOTIFS]
		#self.X_COORD += [n[2]-MOTIFS[0][2]+self.BOOF_GP+self.X_COORD[-1] for n in MOTIFS]
		#self.X_TICKS  = ['P'+str(n) for n in range(len(MOTIFS))]
		#self.X_TICKS += ['G'+str(n) for n in range(len(MOTIFS))]
		#self.X_TICKS  = ['' for n in self.X_TICKS]	# blank them out

		self.PLOT_NUM = {'G':0, 'P':1}

		self.X_COORD  = {'G':[n[2]-MOTIFS[0][2] for n in MOTIFS],
		                 'P':[n[2]-MOTIFS[0][2] for n in MOTIFS]}
		self.X_TICKS  = {'G':['' for n in self.X_COORD['G']],
		                 'P':['' for n in self.X_COORD['P']]}
		self.BOOF_XE  = 100                                   # xlim buffer

		if 'exons' in PLOTTING_META and len(PLOTTING_META['exons']):
			MY_EXONS = PLOTTING_META['exons']
		else:
			MY_EXONS = []
		MY_EXONS = [(n[0]-MOTIFS[0][2], n[1]-MOTIFS[0][2]) for n in MY_EXONS]
		
		if 'xlabel' in PLOTTING_META:
			self.XLAB_STR = PLOTTING_META['xlabel']
		else:
			self.XLAB_STR = {'G':'', 'P':''}

		self.MY_EXNUM = [int((n[0]+n[1])/2) for n in MY_EXONS]
		self.EX_TEXT_OFFSET = (0,70)	# play with these numbers until they look right

		BN = -10000
		BP = 10000
		self.MY_EXON_PATCHES = []
		for n in MY_EXONS:
			polygon = Polygon(np.array([[n[0],BN],[n[0],BP],[n[1],BP],[n[1],BN]]), True)
			self.MY_EXON_PATCHES.append(polygon)

		mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})


	def plot_read(self, pairList, pairList2=None, L=None, rname=None, saveName=None):
		"""
		plot read coordinates vs. observed motif sequence

		pairList  - list of tuples (read_start, read_end, motif_start, motif_end)
		pairList2 - in case you want to plot two reads side by side
		L         - read length
		rname     - read name
		saveName  - name of figure to save
		"""
		if pairList2 == None:
			fig = mpl.figure(0,figsize=(13,5))
			pl  = [pairList]
			spl = [121, 122]
		else:
			fig = mpl.figure(0,figsize=(13,9))
			pl  = [pairList, pairList2]
			spl = [221, 222, 223, 224]

		for i in range(len(pl)):
			plotDat = [[], []]
			for n in pl[i]:
				myLetter = n[2].split('_')[0]
				GP   = self.PLOT_NUM[myLetter]
				inds = int(n[2].split('_')[1])
				indf = int(n[3].split('_')[1])
				#x = [GP*self.L_MOTIF+inds,GP*self.L_MOTIF+indf]
				#x = [self.X_COORD[x[0]], self.X_COORD[x[1]]]
				#y = [n[0],n[1]]
				plotDat[GP].append(([self.X_COORD[myLetter][inds], self.X_COORD[myLetter][indf]], [n[0], n[1]]))
			
			#
			#
			#
			mpl.subplot(spl[i*2])
			for n in plotDat[0]:
				mpl.plot(n[0],n[1],linewidth=4,color='k')

			mpl.xlim(self.X_COORD['G'][0]-self.BOOF_XE,self.X_COORD['G'][-1]+self.BOOF_XE)
			mpl.xticks(self.X_COORD['G'],self.X_TICKS['G'],rotation=90,size=8)
			if L != None:
				mpl.ylim([0,L])
			mpl.grid(linestyle='--')
			mpl.ylabel('read position',fontweight='bold')
			mpl.xlabel(self.XLAB_STR['G'],fontweight='bold')
			if i == 0 and rname != None:
				mpl.suptitle(rname)
			ax = mpl.gca()
			ax.add_collection(PatchCollection(self.MY_EXON_PATCHES, alpha=0.20, color='black'))
			ax.tick_params(axis='x', width=1.5, length=10)
			prevY = ax.get_yticks().tolist()

			#
			#
			#
			mpl.subplot(spl[i*2+1])#, sharey=ax)
			for n in plotDat[1]:
				mpl.plot(n[0],n[1],linewidth=4,color='k')
			
			mpl.xlim(self.X_COORD['P'][0]-self.BOOF_XE,self.X_COORD['P'][-1]+self.BOOF_XE)
			mpl.xticks(self.X_COORD['P'],self.X_TICKS['P'],rotation=90,size=8)
			mpl.yticks(prevY, ['' for n in prevY])
			if L != None:
				mpl.ylim([0,L])
			mpl.grid(linestyle='--')
			mpl.xlabel(self.XLAB_STR['P'],fontweight='bold')
			ax = mpl.gca()
			ax.tick_params(axis='x', width=1.5, length=10)
			ax.tick_params(axis='y', width=0.0, length=0)

		mpl.tight_layout()
		mpl.subplots_adjust(wspace=0.09, top=0.92)
		if saveName != None:
			mpl.savefig(saveName)
		mpl.close(fig)


	def plot_multiple(self, readList, myAlpha=0.0025, maxReads=1000, plotTitle=None, saveName=None):
		"""
		plot read coordinates vs. observed motif sequence, for a list of reads

		readList  - list of: [pairList, readLen, readName]
		myAlpha   - alpha value for each plotted read
		maxReads  - max number of reads to plot (to save time)
		plotTitle - self explanatory
		saveName  - name of figure to save
		"""
		readList = readList[:maxReads]
		fig = mpl.figure(0,figsize=(13,5))
		plotDat = [[], []]
		for i in range(len(readList)):
			[pairList,L,my_rname] = readList[i]
			for n in pairList:
				myLetter = n[2].split('_')[0]
				GP   = self.PLOT_NUM[myLetter]
				inds = int(n[2].split('_')[1])
				indf = int(n[3].split('_')[1])
				#x = [GP*self.L_MOTIF+inds,GP*self.L_MOTIF+indf]
				#x = [self.X_COORD[x[0]], self.X_COORD[x[1]]]
				#y = [n[0],n[1]]
				#plotDat[GP].append(([self.X_COORD[inds], self.X_COORD[indf]], [n[0], n[1]]))
				plotDat[GP].append(([self.X_COORD[myLetter][inds], self.X_COORD[myLetter][indf]], [n[0], n[1]]))
		
		plotting_alpha = max([myAlpha,1./float(len(readList))])

		#mpl.xlim(self.X_COORD[0]-self.BOOF_XE,self.X_COORD[-1]+self.BOOF_XE)
		#mpl.ylim([0,4000])
		#mpl.xticks(self.X_COORD,self.X_TICKS,rotation=90,size=8)
		##mpl.grid(linestyle='--')
		#ax = mpl.gca()
		#ax.add_collection(PatchCollection(self.MY_EXON_PATCHES, alpha=0.20, color='black'))
		#for i in range(len(self.MY_EXNUM)):
		#	mpl.text(self.MY_EXNUM[i]+self.EX_TEXT_OFFSET[0],self.EX_TEXT_OFFSET[1],str(i+1),ha='center')
		#mpl.ylabel('read position',fontweight='bold',fontsize=24)
		#mpl.xlabel(self.XLAB_STR,fontweight='bold',fontsize=24)

		#
		#
		#
		mpl.subplot(121)
		for n in plotDat[0]:
			mpl.plot(n[0],n[1],linewidth=4,color=(0.,0.,0.,plotting_alpha))

		mpl.xlim(self.X_COORD['G'][0]-self.BOOF_XE,self.X_COORD['G'][-1]+self.BOOF_XE)
		mpl.xticks(self.X_COORD['G'],self.X_TICKS['G'],rotation=90,size=8)
		mpl.ylim(bottom=0)
		mpl.grid(linestyle='--')
		mpl.ylabel('read position',fontweight='bold')
		mpl.xlabel(self.XLAB_STR['G'],fontweight='bold')
		if plotTitle != None:
			mpl.suptitle(plotTitle, fontsize=28)
		ax = mpl.gca()
		ax.add_collection(PatchCollection(self.MY_EXON_PATCHES, alpha=0.20, color='black'))
		ax.tick_params(axis='x', width=1.5, length=10)
		prevY = ax.get_yticks().tolist()

		#
		#
		#
		mpl.subplot(122)#, sharey=ax)
		for n in plotDat[1]:
			mpl.plot(n[0],n[1],linewidth=4,color=(0.,0.,0.,plotting_alpha))
		
		mpl.xlim(self.X_COORD['P'][0]-self.BOOF_XE,self.X_COORD['P'][-1]+self.BOOF_XE)
		mpl.xticks(self.X_COORD['P'],self.X_TICKS['P'],rotation=90,size=8)
		mpl.yticks(prevY, ['' for n in prevY])
		mpl.ylim(bottom=0)
		mpl.grid(linestyle='--')
		mpl.xlabel(self.XLAB_STR['P'],fontweight='bold')
		ax = mpl.gca()
		ax.tick_params(axis='x', width=1.5, length=10)
		ax.tick_params(axis='y', width=0.0, length=0)

		mpl.tight_layout()
		mpl.subplots_adjust(wspace=0.09, top=0.89)
		if saveName != None:
			mpl.savefig(saveName)
		mpl.close(fig)


	def plot_AD_hist(self, readList, readList_weak=None, plotTitle=None, saveName=None):
		"""
		FOR DEBUGGING: plot a histogram of gene/pseudogene motif coverage

		readList      - list of: [pairList, readLen, readName]
		readList_weak - same as above, but for low-confidence clusters
		plotTitle - self explanatory
		saveName  - name of figure to save
		"""
		return None
		counts = {}
		for n in readList:
			GP = 1*(n[0][0][2].split('_')[0] == 'G')
			inds = int(n[0][0][2].split('_')[1])
			indf = int(n[0][0][3].split('_')[1])
			(inds, indf) = (min([inds,indf]), max([inds, indf])+1)
			for i in range(inds, indf):
				myT = (GP,i)
				if myT not in counts:
					counts[myT] = 0
				counts[myT] += 1
		x = []
		y = []
		for k in sorted(counts.keys()):
			x.append(self.X_COORD[k[0]*self.L_MOTIF+k[1]])
			y.append(counts[k])

		if readList_weak != None:
			counts2 = {}
			for n in readList_weak:
				GP = 1*(n[0][0][2].split('_')[0] == 'G')
				inds = int(n[0][0][2].split('_')[1])
				indf = int(n[0][0][3].split('_')[1])
				(inds, indf) = (min([inds,indf]), max([inds, indf])+1)
				for i in range(inds, indf):
					myT = (GP,i)
					if myT not in counts2:
						counts2[myT] = 0
					counts2[myT] += 1
			x2 = []
			y2 = []
			for k in sorted(counts2.keys()):
				x2.append(self.X_COORD[k[0]*self.L_MOTIF+k[1]])
				y2.append(counts2[k])

		fig = mpl.figure(0,figsize=(12,6))
		mpl.plot(x,y,linewidth=3)
		if readList_weak != None:
			mpl.plot(x2,y2,linewidth=2,linestyle='--',color='red')
			mpl.legend(['high-conf','weak'],loc=4)
		mpl.xlim(self.X_COORD[0]-self.BOOF_XE,self.X_COORD[-1]+self.BOOF_XE)
		mpl.ylim(bottom=0)
		mpl.xticks(self.X_COORD,self.X_TICKS,rotation=90,size=6)
		mpl.grid(linestyle='--')
		ax = mpl.gca()
		ax.add_collection(PatchCollection(self.MY_EXON_PATCHES, alpha=0.15))
		mpl.ylabel('coverage',fontweight='bold')
		mpl.xlabel(self.XLAB_STR,fontweight='bold')
		if plotTitle != None:
			mpl.title(plotTitle, fontsize=28)
		if saveName != None:
			mpl.savefig(saveName)
		mpl.close(fig)


if __name__ == '__main__':
	pass
