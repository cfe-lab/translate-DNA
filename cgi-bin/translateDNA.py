#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# Random comment
import cgi

form = cgi.FieldStorage()
userinput = form.getvalue("userinput")
runtranslate = form.getvalue("runtranslate")

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
				'D':'ATG', 'B':'TGC', 'N':'ATGC', '-':'ATGC'}

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
def translateDNA(sequence, flag=2):
	sequence = sequence.translate(None, ' \n\r\n').upper()
	aaseq = []
	# Check that the sequence can be divided into codons
	if (len(sequence) % 3 != 0):
		return "Error: sequence length not divisible by 3. Check that you have entire sequence entered"
	# If user wants to output all synonymous mixtures as are, and all non-synonymous as "X"
	if (flag == 2):
		i = 0
		while i < len(sequence):
			codon = resolveCodon(sequence[i:i+3])
			# If the codon has no mixture bases just add it to the amino acid chain
			if len(codon) <= 1:
				aaseq.append(codon_dict[codon[0]])
			# Codon contains mixture base
			else:
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if len(unique) > 1:
					aaseq.append("X")
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
	<body><div class="container">"""
	print "'{}'".format(translateDNA(userinput))
	print "</div></body></html>"
