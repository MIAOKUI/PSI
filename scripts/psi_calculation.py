#!/data/public/tools/python27/bin/python2.7
import sys, itertools, optparse, warnings, re
optParser = optparse.OptionParser(
        usage = "python prepare_exclusion_count.py [options] <flattened_gff_file> <inclusion.count> <exclusion.count> <output_basename>",
        description = "Adding something here",
        epilog = "Adding something here" )

optParser.add_option( "-l", "--readLen", type="int",dest = "readLen",
        default = 100,
        help = "'int number, Indicates the length of short read  (default: 100)" )

optParser.add_option( "-t", "--threshold", type="int",dest = "threshold",
        default = 0,
        help = "'int number, Indicates the threshold of inclusion coverage of an exonic part (default: 0)" )

optParser.add_option( "-v", "--verbose", dest = 'verbose', 
	choices = ( "yes", "no" ),default = 'no', 
	help = "Indicats if print out detail information of annotation (default: yes)) ")

# fetch arguments and options 
(opts, args) = optParser.parse_args()

if len( sys.argv ) == 1:
        optParser.print_help()
        sys.exit(1) 
    
if len( args ) != 4:
        sys.stderr.write( sys.argv[0] + ": Error: Please provide three arguments.\n" )
        sys.stderr.write( "  Call with '-h' to get usage information.\n" )
        sys.exit( 1 ) 
try:
        import HTSeq
except ImportError:
        sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
        sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
        sys.exit(1)

# define exonic part class to store all the related attribute 
class exonicPart: 
	def __init__(self, exonicID):
		self.exonicID = exonicID  	# exonic part id 
		self.geneID = None	 	# gene id 
		self.strand = None		# strand tag +/- 
		self.chr = None			# chromosome number 
		self.start = None		# start coordinate of exonic part  
		self.stop = None		# stop coordinate of exonic part 
		self.exonicLen = None		# The length of exonic part 
		self.inclusion = None 		# inclusion count 
		self.norInclusion = None	# normalized inclusion count 
		self.exclusion = None		# exclusion count 
		self.norExclusion = None	# mormalized exclusion count 
		self.readLen = None		# short read length 
		self.psi = None			# psi 
		self.lowReadHolder = None	# low inclusion read threshold 
		self.lowReadTag = None		# low inclusion read score 
	def getGeneID(self,geneID):
		self.geneID = geneID 
	def getStrand(self, strandTag):
		self.strand = strandTag 
	def getChr(self, chr):
		self.chr = chr 
	def getStart(self, start):
		self.start = start
	def getStop(self, stop): 
		self.stop = stop 
	def getLen(self, length):
		self.exonicLen = length 
	def getInclusion(self, inclusion):
		self.inclusion = inclusion 
	def getExclusion(self, exclusion):
		self.exclusion = exclusion
	def getReadLength(self, readLength):
		self.readLen = readLength
	def getLowReadHolder(self, Threshold):
		self.lowReadHolder = Threshold
		if self.inclusion + self.exclusion < self.lowReadHolder: 
			self.lowReadTag = 0.2
		else: 
			self.lowReadTag = 1 
	# Method calculating psi 
	def getPsi(self):
		self.norInclusion = float(self.inclusion) /(self.readLen-1 + self.exonicLen)
		self.norExclusion = float(self.exclusion) /(self.readLen-1)
		if self.norInclusion + self.norExclusion != 0:  
			self.psi= float(self.norInclusion) /(self.norInclusion + self.norExclusion)
		else: 
			self.psi = 'NA'
	def __repr__(self):
		return "<exonicPart object:" + self.exonicID+">"
	def __str__(self):
		return( "Exonic Part ID:	" + self.exonicID + "\n" 
			 "Exonic Part Start:	" + str(self.start) + "\n"
			 "Exonic Part Stop:	" + str(self.stop) + "\n" 
			 "Inclusion Count:	" + str(self.inclusion) + "\n" 
			 "Exclusion Count:	" + str(self.exclusion) + "\n" 
			 "PSI:			" + str(self.psi) + "\n" ) 
## gtf file parser 
def gtfParser(anotatedLine):
	'''
	Parsing the gtf annotation record and create a list including all annotated information.
	'''
	if anotatedLine.type == 'exonic_part':
		return( [ anotatedLine.attr['gene_id'] + ':' + anotatedLine.attr['exonic_part_number'],			# exonic part id 
			  anotatedLine.attr['gene_id'],	# gene id 
			  anotatedLine.iv.end - anotatedLine.iv.start, # exonic part length 
			  anotatedLine.iv.chrom, 	# chromosome number 
			  anotatedLine.iv.start, 	# exonic part start coordinate 
			  anotatedLine.iv.end,		# exonic part stop coordinate 
			  anotatedLine.iv.strand] )	# strand tag 

# inclusion count file parser
def inclusionParser(inclusionLine):
	'''
	Parsing the inclusion count result return dict object 
	'''
	tmp = re.split(":|\t",inclusionLine.strip('\n'))
	return({'geneID':	tmp[0],
		'exonicID':	tmp[0] + ":" + tmp[1],
		'inclusion':	int(tmp[2]) })

# exclusion count file parser 
def exclusionParser(exclusionLine):
	'''
	Parsing the exclusion  count result return dict object 
	'''
	tmp = re.split(":|\t", exclusionLine.strip('\n'))
	return({'geneID':       tmp[0],
                'exonicID':     tmp[0] + ":" + tmp[1], 
                'exclusion':    int(tmp[2]) }) 

# fetch all the argument 
gtfPath = args[0]
inclusionPath = args[1]
exclusionPath = args[2]
readLength = opts.readLen
threshold = opts.threshold 
basename = args[3]

gtfFile = HTSeq.GFF_Reader(gtfPath,end_included = True)

# read in gtf anotation file 
gtfDict = {}
for gtfline in gtfFile: 
	parsedLine = gtfParser(gtfline)
	if parsedLine is  None:
		continue 
	else: 
		gtfDict[ parsedLine[0] ]= [parsedLine[0],  # exonic ID  
					   parsedLine[1],  # geneID
					   parsedLine[2],  # exonic length 
					   parsedLine[3],  # chrom number 
					   parsedLine[4],  # start 
					   parsedLine[5],  # end 
					   parsedLine[6] ] # strand 

# Read in inclusion and exclusion file and create exonci part object
inclusionFile = open(inclusionPath, 'r')
exclusionFile = open(exclusionPath, 'r')

# iterate the inclusion and exclusion file,calculation psi and output the result 
fout = open(basename+'.psi', 'wt')
for inc,exc in zip(inclusionFile,exclusionFile): 
	if inc[0] == '_':
		break
	 
	incParsed = inclusionParser(inc)
	excParsed = exclusionParser(exc)
	if incParsed['exonicID'] != excParsed['exonicID']:
		print 'Sorting Error: inclusion count and exclusion count may not in a same order !'
		exit(111)
	exonic = exonicPart(incParsed['exonicID'])
	exonic.getGeneID(incParsed[ 'geneID' ]) 
	exonic.getStrand(gtfDict[ exonic.exonicID ][6])
	exonic.getChr(gtfDict[ exonic.exonicID ][3])
	exonic.getStart(gtfDict[ exonic.exonicID ][4])
	exonic.getStop(gtfDict[ exonic.exonicID ][5])
	exonic.getLen(gtfDict[ exonic.exonicID ][2])
	exonic.getInclusion(incParsed[ 'inclusion' ])
	exonic.getExclusion(excParsed[ 'exclusion' ])
	exonic.getReadLength(readLength)
	exonic.getLowReadHolder(threshold)
	exonic.getPsi()
	if opts.verbose == 'yes': 
   		fout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (exonic.exonicID,
									  exonic.psi,
									  exonic.lowReadTag,
									  exonic.inclusion, 
									  exonic.exclusion,
									  exonic.chr, 
									  exonic.start, 
									  exonic.stop, 
									  exonic.strand, 
									  exonic.exonicLen))
	 
   	fout.write( '%s\t%s\t%s\n' % (exonic.exonicID,
				      exonic.psi,
				      exonic.lowReadTag))
fout.close()
