# Checked for Python 3.7

import cgi, re, sys, os

def run(userinput, resolvecharacter, flagselection, highlight, readingframe, button):
	
	# initialization	
	out_str = ""


	##### Function Definitions
	
	
	def parseFasta(lines):
		result = re.findall(r'(>.+[\n\r]+)([^>]+)',lines,re.IGNORECASE)
		tbl = str.maketrans(dict.fromkeys('\n\r\n '))
		result = [y.translate(tbl) for x in result for y in x]
		return result
	
	
	def testWrite():
		currentdir = os.getcwd()
		with open(currentdir+'/tmp/test.txt','w') as f:
			f.write(currentdir)
	
	def printErrors(errors):
		output = ""
		dontrun = False
		if (any(errors.values())):
			output += ('<div class="container word-wrap top-50"><h2>Warning/Error Messages</h2>')
			for e in errors['Error']:
				e = e.split(':')
				if (len(e) > 2) and (e[2] == 'I'):
					output += ("<b>Error: line {}</b>, Sequence '{}' at this line contains illegal characters, translation stopped<br>".format(e[0],e[1]))
					dontrun = True
				else:
					output += ("<b>Error: line {}</b>, Sequence ending at this line is not divisible by 3; last codon has been omitted.<br>".format(e[0]))
			for w in errors['Warning']:
				output += ("Warning: line {} codon {} contains one or two gap characters<br>".format(w.split(':')[0], int(float(w.split(':')[1])))) # +1 to start at 1.
			output += ('</div>')
			return dontrun, output
		elif (dontrun == True):
			return True, output
		else:
			return False, output  


	def checkFasta(fasta, readingframe):
		fasta = fasta.strip()
		errors = {'Warning':[],'Error':[]}
		fasta = fasta.split('\n')
		tbl = str.maketrans(dict.fromkeys('\r '))
		fasta = [x.translate(tbl) for x in fasta]

		# Stack Exchange says this nested function is good!
		def check_errors(seq, index, errors):
			if (seq):
				if re.search('[^ACTGRYKMSWBHDVNX\\-]', seq, re.IGNORECASE):
					errors['Error'].append(str(index)+':'+seq+':I')
			
			# Check the compiled sequence length
			remainder = len(seq) % 3
			if (len(seq[readingframe-1:]) % 3 != 0):
				errors['Error'].append(str(index)+':L')

			problems = checkGaps(seq[readingframe-1:])
			errors['Warning'] = errors['Warning'] + [str(index)+':'+str(x) for x in problems]

		index = 0
		seq = ''
		for line in fasta:
			if len(line) == 0:
				continue
			elif line[0] == '>':
				check_errors(seq, index, errors)
				seq = ''
				index+=2
			else:
				seq += line.upper()
		check_errors(seq, index, errors)  # Check last sequence.
		return errors

	# TODO: make this check for gap characters...	
	def checkSeqs(seqs, readingframe):
		lcount = 1
		errors = {'Warning':[],'Error':[]}
		for seq in seqs:
			if re.search('[^ACTGRYKMSWBHDVNX\\-]', seq, re.IGNORECASE):
				errors['Error'].append(str(lcount)+':'+seq+':I')
			# Check the compiled sequence length
			remainder = len(seq) % 3
			if (len(seq[readingframe-1:]) % 3 != 0):
				errors['Error'].append(str(lcount)+':L')
			
			problems = checkGaps(seq[readingframe-1:])
			errors['Warning'] = errors['Warning'] + [str(lcount)+':'+str(x) for x in problems]

			lcount += 1
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
			  '---':'-', 'XXX':'-', '???':'?'}

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
			return ['???']
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
		#sequence = sequence.translate(None, ' \n\r\n').upper()
		tbl = str.maketrans(dict.fromkeys(' \n\r\n'))  # I think the space is important.
		sequence = sequence.translate(tbl).upper()	
		aaseq = []
		# Check that the sequence can be divided into codons
		#if (len(sequence) % 3 != 0):
		#	return "ERROR"
		# If user wants to output all synonymous mixtures as are, and all non-synonymous as "resolvecharacter"
		i = 0
		while i < len(sequence):
			codon = resolveCodon(sequence[i:i+3])
			if len(codon[0]) < 3: # IF seq is not div by 3 then the last codon will be filtered out with this.
				break # Not continue?

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
	
	
	##### Run analysis

	
	if button == "run":
		out_str += """{% load static %}<!DOCTYPE html><html><head> <title>DNA Translator - Results</title>
			      <link rel="stylesheet" href="{% static "dna_css/style.css" %}">
			      <script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
			      <script src="{% static "dna_script.js" %}"></script> </head><body>"""

		# Accomodate fasta format
		if (userinput[0] == '>'):
			errors = checkFasta(userinput,readingframe)
			is_exit, err_str = printErrors(errors)
			out_str += err_str
			if is_exit:
				return out_str + "</body></html>"  #Is this right?
				
			out_str += ('<div class="spacer"></div><div class="container word-wrap wrapper grey"><div id="sequencestag" unselectable="on" class="noselect">Select all</div><div id="sequences">')
			lines = parseFasta(userinput)
			header = True
			linecount = 1
			out_str += "<table>"
			for line in lines:
				if (header is True):
					out_str += ("<tr id=line{}><td unselectable=\"on\" class=\"noselect linenumber\">{}</td><td>{}</td></tr>".format(linecount, linecount, line))
					header = False
				else:
					remainder = len(line[readingframe-1:]) % 3
					if (remainder == 0):
						out_str += ("<tr id=line{}><td unselectable=\"on\" class=\"noselect linenumber\">{}</td><td>{}</td></tr>".format(linecount, linecount, ('').join(translateDNA(line[readingframe-1:],resolvecharacter,flagselection,highlight))))
					else:
						out_str += ("<tr id=line{}><td unselectable=\"on\" class=\"noselect linenumber\">{}</td><td>{}</td></tr>".format(linecount, linecount, ('').join(translateDNA(line[readingframe-1:-remainder],resolvecharacter,flagselection,highlight))))
					header = True
				linecount += 1
			return out_str + "</table> </div></div>\n</body></html>"

		else:
			tbl = str.maketrans(dict.fromkeys('\r'))
			lines = userinput.translate(tbl).split('\n')
			#lines = userinput.translate(None,'\r').split('\n')
			errors = checkSeqs(lines, readingframe)
			is_exit, err_str = printErrors(errors)
			out_str += err_str
			if is_exit:
				return out_str + "</body></html>"
			
			out_str += ('<div class="spacer"></div><div class="container word-wrap wrapper grey"><div id="sequencestag" unselectable="on" class="noselect">Select all</div><div id="sequences">')
			# Single sequence
			if (len(lines) == 1):
				aa = translateDNA(lines[0][readingframe-1:],resolvecharacter,flagselection,highlight)
				out_str += "<table><td unselectable=\"on\" class=\"noselect linenumber\">1</td>"
				out_str += ("<td>{}</td>".format(('').join(aa)))
				out_str += "</table>"
			# Multiple sequences per line
			elif (len(lines) > 1):
				out_str += ("<table>")
				line_number = 1
				for sequence in lines:
					aa = translateDNA(sequence[readingframe-1:],resolvecharacter,flagselection,highlight)
					out_str += ("<tr><td unselectable=\"on\" class=\"noselect linenumber\">{}</td><td>{}</td></tr>".format(line_number, ('').join(aa)))
					line_number+=1
				out_str += ("</table>")
			return out_str + "</div></div></body></html>"
	else:
		return "Click the translate button to run this analysis.  (Why did I even write a message for this case?)"
