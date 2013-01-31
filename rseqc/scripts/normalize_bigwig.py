#!/usr/bin/env python
'''Normalize bigwig signal to fixed wigsum (equivelent to total reads). Output wiggle file'''
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import string
from optparse import OptionParser

from bx.bbi.bigwig_file import BigWigFile
from qcmodule import BED
from qcmodule import twoList
import numpy
__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.3.3"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

def load_chromsize(file):
	'''read chrom.size file'''
	chromSize={}
	for line in open(file,'r'):
		if line.startswith('#'):continue
		if not line.strip():continue
		fields = line.strip().split()
		chromSize[fields[0]] = int(fields[1])
	return chromSize

def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	
	parser.add_option("-i","--bwfile",action="store",type="string",dest="BigWig_File",help="Input BigWig file. [required]")
	parser.add_option("-o","--output",action="store",type="string",dest="output_wig",help="Output wig file. [required]")
	parser.add_option("-s","--chromSize",action="store",type="string",dest="chromSize",help="Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome. [required]")
	parser.add_option("-t","--wigsum",action="store",type="int",dest="total_wigsum",default=100000000,help="Specified wigsum. 100000000 equals to coverage of 1 million 100nt reads. default=%default  [optional]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed format. [optional]")	
	parser.add_option("-c","--chunk",action="store",type="int",dest="chunk_size",default=100000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp) [optional]")
	(options,args)=parser.parse_args()
	
	if not (options.BigWig_File and options.output_wig and options.chromSize):
		parser.print_help()
		sys.exit(0)

	OUT=open(options.output_wig,'w')
	bw = BigWigFile( file=open(options.BigWig_File) )
	chrom_sizes = load_chromsize(options.chromSize)	
	exons=[]
	WIG_SUM=0.0
	if (options.refgene_bed):	
		print >>sys.stderr, "Extract exons from " + options.refgene_bed
		obj = BED.ParseBED(options.refgene_bed)
		exons = obj.getExon()
		print >>sys.stderr, "Merge overlapping exons ..."
		exons = BED.unionBed3(exons)
		print >>sys.stderr, "Calculate wigsum covered by " + options.refgene_bed + ' only'
		for chrom,st,end in exons:
			try: bw.get_as_array(chrom,0,1).size
			except:continue

			bw_signal = bw.get_as_array(chrom,st,end)
			tmp = numpy.nansum(bw_signal)			#nan will be ignored. but if all items are 'nan', the result summay is 'nan' NOT 0
			if numpy.isnan(tmp):continue	
			WIG_SUM += tmp
		print >>sys.stderr, "Total wigsum is %.2f\n" % WIG_SUM
	else:
		print >>sys.stderr, "Calculate wigsum from " + options.BigWig_File
		for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom

			try: bw.get_as_array(chr_name,0,1).size
			except:
				print >>sys.stderr, "Skip " + chr_name + "!"
				continue

			print >>sys.stderr, "Processing " + chr_name + " ..."	
			for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
				bw_signal = bw.get_as_array(interval[0],interval[1],interval[2])
				tmp = numpy.nansum(bw_signal)
				if numpy.isnan(tmp):continue
				WIG_SUM += tmp
		print >>sys.stderr, "\nTotal wigsum is %.2f\n" % WIG_SUM
	
	try:
		weight = options.total_wigsum/WIG_SUM
	except:
		"Error, WIG_SUM cannot be 0"
		eys.exit(1)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	print >>sys.stderr, "Normalizing bigwig file, output wiggle file"
	for chr_name, chr_size in chrom_sizes.items():          #iterate each chrom
		

		try: bw.get_as_array(chr_name,0,1).size
		except:
			print >>sys.stderr, "Skip " + chr_name + "!"
			continue
		
		print >>sys.stderr, "Writing " + chr_name + " ..."
		OUT.write('variableStep chrom='+chr_name+'\n')
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
			coord = interval[1]
			bw_signal = bw.get_as_array(chr_name,interval[1],interval[2])
			tmp = numpy.nansum(bw_signal)
			if numpy.isnan(tmp):continue
			bw_signal = numpy.nan_to_num(bw_signal)
			for v in bw_signal:
				coord +=1
				if v != 0: print >>OUT, "%d\t%.4f" % (coord,v*weight)
			
			
if __name__=='__main__':
	main()
