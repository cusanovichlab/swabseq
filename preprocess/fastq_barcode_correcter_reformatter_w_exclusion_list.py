#import pysam #v0.8.1
import sys
import gzip
import io
import Levenshtein
import subprocess
import itertools
import more_itertools

infastq1 = sys.argv[1]
infastq2 = sys.argv[2]
infastq3 = sys.argv[3]
inbarcodes = sys.argv[4]
inminuslist = sys.argv[5]
distmetric = sys.argv[6]
outprefix = sys.argv[7]

print("Fastq 1 input: " + infastq1)
print("Fastq 2 input: " + infastq2)
print("Fastq 3 input: " + infastq3)
print("Samplesheet input: " + inbarcodes)
print("minuslist input: " + inminuslist)
print("Distance metric: " + distmetric)
print("Output prefix: " + outprefix)

#This version of editcheck just does a light cleaning - looking to see if ANY reference barcodes are within 1 edit of the sequenced barcode.
#def editcheck(barc,reflist):
#	winner_ed = 10
#	for barcode in reflist.keys():
#		curred = Levenshtein.distance(barc,barcode)
#		if curred < winner_ed:
#			winner_ed = curred
#			if winner_ed == 1:
#				break
#	return(winner_ed)

#This version of editcheck does a much more selective search, checking to see if there are ties for the best matching barcode and returning the sequence.
def editcheck(barc,reflist):
	winner_ed = 10
	runnerup_ed = 10
	for barcode in reflist.keys():
		curred = Levenshtein.distance(barc,barcode)
		if curred <= winner_ed:
			runnerup_ed = winner_ed
			winner = barcode
			winner_ed = curred
	#if winner_ed > 1:
	#	winner = '_CTF' + '_'*(len(barc)-4)
	if runnerup_ed - winner_ed < 1:
		winner_ed = 10
	#	winner = '_AMBIG' + '_'*(len(barc)-6)
	return(winner,winner_ed)

#This code modified from Andrew Hill sciATAC pipeline
def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of pluslists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for n_mismatch,mismatch_pluslist in enumerate(mismatch_map):
        corrected = mismatch_pluslist.get(barcode, None)

        if corrected:
            return(corrected,n_mismatch)

    return(None,10)


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    Yield:
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_pluslist_map(pluslist, minuslist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode pluslist.
    Args:
        pluslist (set of str): set of pluslist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their pluslist sequences
    """

    mismatch_to_pluslist_map = [None] * (edit_distance + 1)

    mismatch_to_pluslist_map[0] = {k: k for k in pluslist}

    if minuslist == "None":
    	plusminuslist = pluslist
    else:
    	plusminuslist = pluslist + minuslist

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(plusminuslist)

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_pluslist_map[mismatch_count] = {}
        if minuslist != "None":
	        for sequence in minuslist:
	            sequence = sequence.upper()

	        	# Generate all possible mismatches in range
	            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

	            # Construct a mapping to the intended sequences
	            for mismatch in mismatches:
	                # Track all minuslist matches for removal later
	                if mismatch not in mismatch_to_pluslist_map[mismatch_count]:
	                    conflicting_mismatches.append(mismatch)
	                mismatch_to_pluslist_map[mismatch_count][mismatch] = sequence

        for sequence in plusminuslist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_pluslist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_pluslist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_pluslist_map[mismatch_count]:
                del mismatch_to_pluslist_map[mismatch_count][mismatch]

    return mismatch_to_pluslist_map

#End of code modified from Andrew Hill sciATAC pipeline

print("Building barcode list...")
inindex = open(inbarcodes,'r')
plater = 0
indexset = {}
index2set = {}
comboset = {}
for line in inindex:
	liner = line.strip().split(",")
	if len(liner) == 0:
		continue
	if plater == 0:
		if liner[0] == "Plate_ID":
			indexind = liner.index("index")
			index2ind = liner.index("index2")
			plater = 1
			continue
		else:
			continue
	try:
		indexset[liner[indexind]]
	except KeyError:
		indexset[liner[indexind]] = ""
	try:
		index2set[liner[index2ind]]
	except KeyError:
		index2set[liner[index2ind]] = ""
	try:
		comboset[liner[indexind] + liner[index2ind]]
	except KeyError:
		comboset[liner[indexind] + liner[index2ind]] = ""

inindex.close()
print("Found " + str(len(comboset.keys())) + " barcode combinations, including " + str(len(indexset.keys())) + " P7 indices and " + str(len(index2set.keys())) + " P5 indices...")

if inminuslist != "None":
	print("Building minuslist...")
	inminus = open(inminuslist,'r')
	minusindexset = {}
	minusindex2set = {}
	minuscomboset = {}
	#combocount = 0
	for line in inminus:
		liner = line.strip().split("\t")
		try:
			#combocount += 1
			indexset[liner[0]]
		except KeyError:
			try:
				minusindexset[liner[0]]
			except KeyError:
				minusindexset[liner[0]] = ""
		try:
			index2set[liner[1]]
		except KeyError:
			try:
				minusindex2set[liner[1]]
			except KeyError:
				minusindex2set[liner[1]] = ""

		try:
			comboset[liner[0] + liner[1]]
			print("Error! Found conflict between minuslist and pluslist barcodes! Offending barcode: " + liner[0] + liner[1])
			sys.exit()
		except KeyError:
			minuscomboset[liner[0] + liner[1]] = ""

	inindex.close()
	print("Found " + str(len(minuscomboset.keys())) + " minuslist barcode combinations, including " + str(len(minusindexset.keys())) + " P7 indices not in pluslist and " + str(len(minusindex2set.keys())) + " P5 indices not in pluslist...")
else:
	print("No minuslist provided...")

if distmetric == "Hamming":
	print("Building P7 mismatch map...")
	if inminuslist == "None":
		minusindexset =  {"None": ""}
		minusindex2set =  {"None": ""}
	p7_correction_map = construct_mismatch_to_pluslist_map(list(indexset.keys()), list(minusindexset.keys()), 1)
	print("P7 mismatch map links " + str(len(p7_correction_map[1].keys())) + " mismatches to " + str(len(set(p7_correction_map[1].values()))) + " valid bacodes")
	print("Examples: " + str(more_itertools.take(5, p7_correction_map[1].items())))
	print("Building P5 mismatch map...")
	p5_correction_map = construct_mismatch_to_pluslist_map(list(index2set.keys()), list(minusindex2set.keys()), 1)
	print("P5 mismatch map links " + str(len(p5_correction_map[1].keys())) + " mismatches to " + str(len(set(p5_correction_map[1].values()))) + " valid bacodes")
	print("Examples: " + str(more_itertools.take(5, p5_correction_map[1].items())))

#monkeyreads

totreads = 0
exactmatch = 0
editmatch = 0
failed = 0
failedmatch = 0
#Option 1 - doesn't seem to work in python 3
#inner1 = gzip.open(infastq1,"rb")
#inner2 = gzip.open(infastq2,"rb")
#inner3 = gzip.open(infastq3,"rb")
#readsin1 = io.BufferedReader(inner1)
#readsin2 = io.BufferedReader(inner2)
#readsin3 = io.BufferedReader(inner3)
#Option 2 - works, but SLOW...
#readsin1 = gzip.open(infastq1,"r")
#readsin2 = gzip.open(infastq2,"r")
#readsin3 = gzip.open(infastq3,"r")
#Option 3 - wasteful with space, but should be faster...
print("Unpacking read 1...")
subprocess.call(["gunzip",infastq1])
print("Unpacking index 1...")
subprocess.call(["gunzip",infastq2])
print("Unpacking index 2...")
subprocess.call(["gunzip",infastq3])
print("Processing reads...")
if distmetric == "Levenshtein" and inminuslist != "None":
	print("Merging minuslist into pluslist as Levenshtein distance was selected...")
	indexset.update(minusindexset)
	index2set.update(minusindex2set)

readsin1 = open(infastq1.split(".gz")[0],"r")
readsin2 = open(infastq2.split(".gz")[0],"r")
readsin3 = open(infastq3.split(".gz")[0],"r")
outfastq1 = open(outprefix + ".R1.fastq",'w')
outfastq2 = open(outprefix + ".I1.fastq",'w')
outfastq3 = open(outprefix + ".I2.fastq",'w')
for line in readsin1:
	totreads += 1
	if totreads % 5000000 == 0:
		print(str(totreads) + " total reads processed...")
		print(str(editmatch+exactmatch) + " reads successfullly demultiplexed so far...")
	liner = line.strip().split()
	dump1 = next(readsin2).strip()
	#barc1 = next(readsin2).strip().decode()
	barc1 = next(readsin2).strip()
	dump2 = next(readsin3).strip()
	#barc2 = next(readsin3).strip().decode()
	barc2 = next(readsin3).strip()
	#print(barc1 + " + " + barc2)
	try:
		indexset[barc1]
		#b1ed = 0
		b1ed = [barc1,0]
		#print("barc1 edit distance = 0")
	except KeyError:
		if distmetric == "Hamming":
			b1ed = correct_barcode(barc1,p7_correction_map)
		if distmetric == "Levenshtein":
			b1ed = editcheck(barc1,indexset)
		#print("barc1 edit distance = " + str(b1ed[1]))
		#if b1ed > 1:
		if b1ed[1] > 1:
			failed += 1
			next(readsin1).strip()
			next(readsin1).strip()
			next(readsin1).strip()
			next(readsin2).strip()
			next(readsin2).strip()
			next(readsin3).strip()
			next(readsin3).strip()
			continue
	try:
		index2set[barc2]
		#b2ed = 0
		b2ed = [barc2,0]
		#print("barc2 edit distance = 0")
	except KeyError:
		if distmetric == "Hamming":
			b2ed = correct_barcode(barc2,p5_correction_map)
		if distmetric == "Levenshtein":
			b2ed = editcheck(barc2,index2set)
		#print("barc2 edit distance = " + str(b2ed[1]))
		#if b2ed > 1:
		if b2ed[1] > 1:
			failed += 1
			next(readsin1).strip()
			next(readsin1).strip()
			next(readsin1).strip()
			next(readsin2).strip()
			next(readsin2).strip()
			next(readsin3).strip()
			next(readsin3).strip()
			continue
	try:
		comboset[b1ed[0] + b2ed[0]]
	except KeyError:
		failedmatch += 1
		#print(b1ed[0] + b2ed[0])
		#print(b1ed)
		#print(b2ed)
		next(readsin1).strip()
		next(readsin1).strip()
		next(readsin1).strip()
		next(readsin2).strip()
		next(readsin2).strip()
		next(readsin3).strip()
		next(readsin3).strip()
		continue
	#if b1ed == 0 and b2ed == 0:
	if b1ed[1] == 0 and b2ed[1] == 0:
		exactmatch += 1
	#if b1ed == 1 or b2ed == 1:
	if b1ed[1] == 1 or b2ed[1] == 1:
		editmatch += 1
#	outfastq1.write(liner[0].decode() + " 1:N:0:" + barc1 + "+" + barc2 + "\n")
#	outfastq1.write(next(readsin1).strip().decode()[0:26] + "\n")
#	outfastq1.write(next(readsin1).strip().decode() + "\n")
#	outfastq1.write(next(readsin1).strip().decode()[0:26] + "\n")
#	outfastq2.write(liner[0].decode() + " 1:N:0:" + barc1 + "+" + barc2 + "\n")
#	outfastq2.write(barc1 + "\n")
#	outfastq2.write(next(readsin2).strip().decode() + "\n")
#	outfastq2.write(next(readsin2).strip().decode() + "\n")
#	outfastq3.write(liner[0].decode() + " 2:N:0:" + barc1 + "+" + barc2 + "\n")
#	outfastq3.write(barc2 + "\n")
#	outfastq3.write(next(readsin3).strip().decode() + "\n")
#	outfastq3.write(next(readsin3).strip().decode() + "\n")
	outfastq1.write(liner[0] + " 1:N:0:" + barc1 + "+" + barc2 + "\n")
	outfastq1.write(next(readsin1).strip()[0:26] + "\n")
	outfastq1.write(next(readsin1).strip() + "\n")
	outfastq1.write(next(readsin1).strip()[0:26] + "\n")
	outfastq2.write(liner[0] + " 1:N:0:" + barc1 + "+" + barc2 + "\n")
	outfastq2.write(barc1 + "\n")
	outfastq2.write(next(readsin2).strip() + "\n")
	outfastq2.write(next(readsin2).strip() + "\n")
	outfastq3.write(liner[0] + " 2:N:0:" + barc1 + "+" + barc2 + "\n")
	outfastq3.write(barc2 + "\n")
	outfastq3.write(next(readsin3).strip() + "\n")
	outfastq3.write(next(readsin3).strip() + "\n")

outfastq1.close()
outfastq2.close()
outfastq3.close()

print("Total reads processed: " + str(totreads))
print("Exact match barcodes processed: " + str(exactmatch))
print("Edit match barcodes processed: " + str(editmatch))
print("Unmatched barcodes processed: " + str(failed))
print("Invalid barcode pairs processed: " + str(failedmatch))
print("Tiddying up...")

print("Packing up read 1...")
subprocess.call(["gzip","-f",outprefix + ".R1.fastq"])
subprocess.call(["gzip",infastq1.split(".gz")[0]])
print("Packing up index 1...")
subprocess.call(["gzip","-f",outprefix + ".I1.fastq"])
subprocess.call(["gzip",infastq2.split(".gz")[0]])
print("Packing up index 2...")
subprocess.call(["gzip","-f",outprefix + ".I2.fastq"])
subprocess.call(["gzip",infastq3.split(".gz")[0]])
print("Done.")
