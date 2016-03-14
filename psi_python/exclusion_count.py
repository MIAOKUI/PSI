#!/bin/env python
import sys, itertools, optparse, warnings

try:
	import HTSeq
except ImportError:
	sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
	sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )



# define junction parser function
def junction_parser(junctionsLine):
	'''
	This function accept each line of tophat output junction file
	and parsed them into genomic interval object defined in htseq library.
	The genomic interval object can't be used for coverage counting.
	'''
	tmp = junctionsLine.split()
	gapWidth = tmp[10].split(",")
	chr = tmp[0]
	strand = tmp[5]
	start = int(tmp[1])+ int(gapWidth[0])
	stop = int(tmp[2]) - int(gapWidth[1])
	count = int(tmp[4])
	iv = HTSeq.GenomicInterval(chr,start,stop,strand)
	return [iv, count]

# define junction filtering function
def junction_filter(junctionIv, gtfStart, gtfStop):
	'''
	This function will filter the junction record and check if the junction
	range is valid. And return a boolean value ( True, False) to indicate the 
	if the junction record pass the filter or not . The start coordinate 
	of a junction should be at least one end coordinate of exonic part in gtf file. 
	And the end coordinate of a junction should be ad least one start coordinate of 
	exonic part in gtf file.   
	'''
	leftBoundary = (str(junctionIv.chrom) + ':' + str(junctionIv.start)) in gtfStop
	rightBoundary = (str(junctionIv.chrom) + ':' + str(junctionIv.end)) in gtfStart 
	return( leftBoundary and rightBoundary )


def main():
	optParser = optparse.OptionParser(
		usage = "python prepare_exclusion_count.py [options] <flattened_gff_file> <junctions.bed> <output_basename>",
		description = "Adding something here",
		epilog = "Adding something here" )

	optParser.add_option( "-f", "--filtered", type="choice",dest = "filtered",
		choices = ( "no", "yes" ), default = "yes",
		help = "'yes' or 'no'. Indicates whether filtering junction file based on present annotation (default: yes)" )

	if len( sys.argv ) == 1:
		optParser.print_help()
		sys.exit(1) 

	# fetch arguement and option	
	(opts, args) = optParser.parse_args()
	if len( args ) != 3:
		sys.stderr.write( sys.argv[0] + ": Error: Please provide three arguments.\n" )
		sys.stderr.write( "  Call with '-h' to get usage information.\n" )
		sys.exit( 1 )
		
	# Import gtf file
	gtfFile = HTSeq.GFF_Reader(args[0], end_included = True)

	# Open junction file 
	junctionsFile = open(args[1]) 

	# Get out file basename 
	outFile = args[2]



	# Open gtf file 
	features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

	# initial exclusion count dict 
	counts = {}

	# initial gtf exonic part start and end filter checking list
	gtfStart = [] 
	gtfStop = []

	# iterate the whole gtf file, initiate the exclusion count, create filter checking list 
	for f in  gtfFile:
		if f.type == "exonic_part":
			f.name = f.attr['gene_id'] + ":" + f.attr['exonic_part_number']
			features[f.iv] += f
			counts[f.name] = 0
			if opts.filtered == 'yes':
				gtfStart.append(str(f.iv.chrom)+ ':' + str(f.iv.start))
				gtfStop.append(str(f.iv.chrom) + ':' + str(f.iv.end))
	

	# Iterate the junction file and generate exclusion count
	for junctionsLine in junctionsFile: 
		# skip junction.bed file header 
		if junctionsLine[0:5] == "track":
			continue
	 	# parse junction record line 
		tmp = junction_parser(junctionsLine)
	
		# filtering 
		if opts.filtered == 'yes': 
			if not junction_filter(tmp[0],gtfStart, gtfStop): 
				continue 
		# counting exclusion 
		for iv,step_set in features[tmp[0]].steps():
			for hit in list(step_set):
				if tmp[0].start < hit.iv.start and tmp[0].end > hit.iv.end:  
					counts[hit.name] += tmp[1]
	junctionsFile.close() 

	# Write result 
	fout = open(outFile+".exclusion","w")
	for name in sorted(counts.keys()):
		fout.write( "%s\t%d\n" % ( name, counts[name] ) )
	fout.close()

if __name__ == '__main__':
	main()
