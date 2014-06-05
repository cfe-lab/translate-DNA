#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import cgi, re

form = cgi.FieldStorage()
userinput = form.getvalue("userinput")
runtranslate = form.getvalue("runtranslate")
flagselection = int(form.getvalue("flagselection"))
resolvecharacter = form.getvalue("resolvecharacter")
highlight = form.getvalue("highlight")
if (resolvecharacter is None):
	resolvecharacter = 'X'

def parseFasta(lines):
	result = re.findall('(>.+[\\n\\r\\n])([\\*\\-ACTGRYKMSWBHDVN\\:\\n\\r\\n]+)', lines, re.IGNORECASE)
	result = [y.translate(None, '\n\r\n') for x in result for y in x]
	return result

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
	if (len(sequence) % 3 != 0):
		return "Error: sequence length not divisible by 3. Check that you have entire sequence entered"
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
	</head>
	<body><div class="container word-wrap">"""
	#print "flag selection: {}<br>resolve character: {}".format(flagselection,resolvecharacter)
	# To accomodate fasta format
	if (userinput[0] == '>'):
		lines = parseFasta(userinput)
		#print repr(lines)
		header = True
		for line in lines:
			if (header is True):
				print "{}<br>".format(line)
				header = False
			else:
				print "{}<br>".format(('').join(translateDNA(line,resolvecharacter,flagselection,highlight)).upper())
				header = True
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
	print "</div></body></html>"
