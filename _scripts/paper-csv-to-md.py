# convert bibtex files to markdown files for website

# example output from Bedford lab
# layout: paper
# title: Optimization of gene expression by natural selection
# image: /images/papers/bedford-optimization.png
# authors: Bedford T, Hartl DL.
# year: 2009
# ref: Bedford and Hartl. 2009. Proc Natl Acad Sci USA.
# journal: "Proc Natl Acad Sci USA 106: 1133-1138."
# pdf: /pdfs/papers/bedford-optimization.pdf
# supplement: /pdfs/papers/bedford-optimization-supp.pdf
# doi: 10.1073/pnas.0812009106 

# starting with csv output from Zotero
papers = 'BrazPubs.csv'

import pandas as pd
df = pd.read_csv(papers)

for index, row in df.iterrows():
	
	# the new filename must begin with date YYYY-MM-DD to be recognized by website
	# so the "Date" field in the csv file should be formatted accordingly
	date = row["Date"]
	if len(date) < 10: date = date + '-01'
	if len(date) < 10: date = date + '-01'
	if len(date) < 10:
		print("error: incorrect date format")
		break
	
	# want the following fields from the csv:
	title = row["Title"]
	authors = row["Author"]
	first = str(authors.split(',')[0])
	year = str(row["Publication Year"])
	journal = str(row["Publication Title"])
	abbrev = str(row["Journal Abbreviation"])
	doi = str(row["DOI"])
	abstract = str(row["Abstract Note"])
	pdf = row["File Attachments"].split('/')[-1]
	image = pdf.replace('.pdf', '.png').replace('.PDF', '.png')
	url = row["Url"]
	url = '[' + url + '](' + url + ')' 
	
	# create reference from fields above
	ref = first + " et al. " + year + ' ' + abbrev
	
	# create a new filename for each paper:
	sep = '-'
	filename = sep.join([first,abbrev]) + '.md'
	filename = date + '-' + filename
	#pdf = filename.replace('.md', '.pdf')
	
	with open(filename, 'w') as o:
		o.write('---\n')
		o.write('layout: paper\n')
		
		o.write('title: >-\n')	# the >- characters escape the : and - in the title, avoiding YAML errors
		o.write('    ')			# indent with four spaces
		o.write(title)
		o.write('\n')
		
		o.write('image: /images/papers/')
		o.write(image)
		o.write('\n')
				
		o.write('authors: ')
		o.write(authors)
		o.write('\n')
		
		o.write('year: ')
		o.write(year)
		o.write('\n')
		
		o.write('ref: ')
		o.write(ref)
		o.write('\n')		
		
		o.write('journal: >-\n')	# the >- characters escape the : and - in the journal, avoiding YAML errors
		o.write('    ')			# indent with four spaces
		o.write(journal)
		o.write('\n')
		
		o.write('pdf: /pdfs/papers/')
		o.write(pdf)
		o.write('\n')
		
		o.write('supplement: /pdfs/supplements/')
		o.write(pdf)
		o.write('\n')
		
		o.write('doi: ')
		o.write(doi)
		o.write('\n')
		
		o.write('---\n\n')
		
		o.write('Link to online paper: ')
		o.write(url)
		o.write('\n')
		
		o.write('\n# Abstract\n\n')
		o.write(abstract)
		o.write('\n\n')
		