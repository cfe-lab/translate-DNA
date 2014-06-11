#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import cgi, re, sys, os

form = cgi.FieldStorage()
userinput = form.getvalue("userinput")
runtranslate = form.getvalue("runtranslate")
flagselection = int(form.getvalue("flagselection"))
resolvecharacter = form.getvalue("resolvecharacter")
highlight = form.getvalue("highlight")
readingframe = int(form.getvalue("readingframe"))
if (resolvecharacter is None):
	resolvecharacter = 'X'

def parseFasta(lines):
	#result = re.findall('(>.+[\\n\\r\\n])([\\*\\-ACTGRYKMSWBHDVN\\:\\n\\r\\n]+)', lines, re.IGNORECASE)
	result = re.findall(r'(>.+[\n\r]+)([^>]+)',lines,re.IGNORECASE)
	result = [y.translate(None, '\n\r\n') for x in result for y in x]
	return result

def testWrite():
	currentdir = os.getcwd()
	with open(currentdir+'/tmp/test.txt','w') as f:
		f.write(currentdir)

def printErrors(errors):
	if (any(errors.values())):
		print '<div class="container word-wrap top-50"><h2>Warning/Error Messages</h2>'
		for e in errors['Error']:
			print "Error: line {}, Sequence ending at this line is not divisible by 3<br>".format(e.split(':')[0])
		for w in errors['Warning']:
			print "Warning: line {} codon {} contains one or two gap characters<br>".format(w.split(':')[0],w.split(':')[1])
		print '</div>'
	else:
		return False

# Error: 'X:L' means line X has sequence not divisible by 3
# Warning: 'X:Y' means line X codon Y has a gap character
def checkInput(input,readingframe):
	input = input.strip()
	errors = {'Warning':[],'Error':[]}
	# If input is fasta format
	if (input[0] == '>'):
		#input = parseFasta(input)
		input = input.split('\n')
		input = [x.translate(None, '\r ') for x in input]
		seq = ''
		for line in enumerate(input):
			if line[1][0] != '>':
				seq += line[1]
				if (readingframe == 1):
					problems = checkGaps(line[1])
					errors['Warning'] = errors['Warning'] + [str(line[0]+1)+':'+str(x) for x in problems]
				if (readingframe == 2):
					problems = checkGaps(line[1][1:])
					errors['Warning'] = errors['Warning'] + [str(line[0]+1)+':'+str(x) for x in problems]
				if (readingframe == 3):
					problems = checkGaps(line[1][2:])
					errors['Warning'] = errors['Warning'] + [str(line[0]+1)+':'+str(x) for x in problems]
			else:
				#check seq for proper length
				if (len(seq) % 3 != 0):
					errors['Error'].append(str(line[0])+':L')
				seq = ''
		#check last seq for proper length
		if (len(seq) % 3 != 0):
			errors['Error'].append(str(line[0])+':L')
	return errors

def checkFasta(fasta, readingframe):
	fasta = fasta.strip()
	errors = {'Warning':[],'Error':[]}
	fasta = fasta.split('\n')
	fasta = [x.translate(None, '\r ') for x in fasta]
	seq = ''
	for line in enumerate(fasta):
		if (line[1][0] == '>'):
			# Check the compiled sequence length
			remainder = len(line) % 3
			if (len(seq[readingframe-1:-remainder]) % 3 != 0):
				errors['Error'].append(str(line[0])+':L')
			seq = ''
			continue
		seq += line[1]
		problems = checkGaps(line[1][readingframe-1:])
		errors['Warning'] = errors['Warning'] + [str(line[0]+1)+':'+str(x) for x in problems]
	return errors

# 0: sequence is not divisible by 3
# [A,B,C,...]: gap character detected at codon A, B, C, ...
def checkGaps(seq):
	seq = seq.upper()
	i = 0
	gaps = []
	while i < len(seq):
		codon = seq[i:i+3]
		xes = codon.count('X')
		dashes = codon.count('-')
		if (xes + dashes != 3) and ((1 <= xes < 3) or (1 <= dashes < 3)):
			gaps.append(i/3 + 1)
		i += 3
	return gaps

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
			  'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
			  'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
			  'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
			  'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
			  'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
			  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
			  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
			  'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
			  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
			  'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
			  'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
			  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
			  'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
			  'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
			  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
			  '---':'-', 'XXX':'?'}

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT',
				'S':'CG', 'M':'AC', 'V':'AGC', 'H':'ATC',
				'D':'ATG', 'B':'TGC', 'N':'ATGC', '-':'-'}

# This amazing little function will take any codon and return a list of all possible
# codons it could resolve to (only one if there are no mixtures)
def resolveCodon(codon):
	nonmix = []
	if (codon in codon_dict):
		return [codon]
	elif (codon.count('-') + codon.count('X') == 3):
		return ['---']
	elif (1 <= codon.count('-') <= 2) or (1 <= codon.count('X') <= 2):
		return ['XXX']
	for base in codon:
		# Check for mixtures
		if (base in mixture_dict):
			if (not nonmix):
				nonmix = [x for x in mixture_dict[base]]
			else:
				nonmix = [x+y for x in nonmix for y in mixture_dict[base]]
		else:
			if (not nonmix):
				nonmix.append(base)
			else:
				nonmix = [x+base for x in nonmix]
	return nonmix

# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter, flag=2, highlight='True'):
	sequence = sequence.translate(None, ' \n\r\n').upper()
	aaseq = []
	# Check that the sequence can be divided into codons
	#if (len(sequence) % 3 != 0):
	#	return "ERROR"
	# If user wants to output all synonymous mixtures as are, and all non-synonymous as "resolvecharacter"
	i = 0
	while i < len(sequence):
		codon = resolveCodon(sequence[i:i+3])
		# If the codon has no mixture bases just add it to the amino acid chain
		if len(codon) <= 1:
			aaseq.append(codon_dict[codon[0]])
		# Codon contains mixture base
		else:
			# If flag is set to 1
			if (flag == 1):
				if (highlight == 'True'):
					aaseq.append('<span class="highlight">'+resolvecharacter+'</span>')
				else:
					aaseq.append(resolvecharacter)
			# If flag is set to 2
			elif (flag == 2):
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					if (highlight == 'True'):
						aaseq.append('<span class="highlight">'+resolvecharacter+'</span>')
					else:
						aaseq.append(resolvecharacter)
				else:
					aaseq.append(unique.pop())
			# If flag is set to 3
			else:
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					if (highlight == 'True'):
						aaseq.append('<span class="highlight">['+('/').join(unique)+']</span>')
					else:
						aaseq.append('['+('/').join(unique)+']')
				else:
					aaseq.append(unique.pop())
		i += 3
	return aaseq

if (runtranslate is not None):
	print "Content-Type: text/html"
	print
	print """<!DOCTYPE html><html>
	<head>
	<link rel="stylesheet" href="../css/style.css">
	<script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
	<script src="../cgi-bin/script.js"></script>
	</head>
	<body>"""
	#print "flag selection: {}<br>resolve character: {}".format(flagselection,resolvecharacter)
	errors = checkFasta(userinput,readingframe)
	printErrors(errors)
	print '<div class="spacer"></div><div class="container word-wrap wrapper grey"><div id="sequencestag">Select all</div><div id="sequences">'
	# To accomodate fasta format
	if (userinput[0] == '>'):
		lines = parseFasta(userinput)
		#print repr(lines)
		header = True
		linecount = 1
		for line in lines:
			if (header is True):
				print "<span id=line{}>{}</span><br>".format(linecount,line)
				header = False
			else:
				# Reading frame 1
				if (readingframe == 1):
					remainder = len(line) % 3
				# Reading frame 2
				elif (readingframe == 2):
					remainder = len(line[1:]) % 3
				# Reading frame 3
				else:
					remainder = len(line[2:]) % 3
				if (remainder == 0):
					print "<span id=line{}>{}</span><br>".format(linecount,('').join(translateDNA(line[readingframe-1:],resolvecharacter,flagselection,highlight)))
				else:
					print "<span id=line{}>{}</span><br>".format(linecount,('').join(translateDNA(line[readingframe-1:-remainder],resolvecharacter,flagselection,highlight)))
				header = True
			linecount += 1
		sys.exit()
	lines = userinput.translate(None,'\r').split('\n')
	# Single sequence
	if (len(lines) == 1):
		aa = translateDNA(lines[0],resolvecharacter,flagselection,highlight)
		print "{}".format(('').join(aa))
	# Multiple sequences per line
	elif (len(lines) > 1):
		print "<table>"
		for sequence in lines:
			aa = translateDNA(sequence,resolvecharacter,flagselection,highlight)
			print "<tr><td>{}</td></tr>".format(('').join(aa))
		print "</table>"
	print "</div></div></body></html>"
