# usage: python decipher.py input.txt output.txt
import urllib
import urllib2
import sys,re
from lxml import html
import os

def parse_ucsc(link):
	"""parse ucsc link for decipher"""
	pattern = re.compile(".+HREF='.+&c=(\w+)&o=(\d+)&t=(\d+).+&i=(\d+)'.+TITLE='(.+)'")
	item = re.match(pattern,link)
	pid = item.group(4)
	pos = item.group(1)+':'+item.group(2)+'-'+item.group(3)
	phenotype = item.group(5)
	if pid == phenotype:
		phenotype = '-'
	return [pid,pos,phenotype]

ucsc_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks'
values = {'org' : 'human',
		  'db': 'hg19',
		  'decipher': 'full'}

decipher_url = 'https://decipher.sanger.ac.uk/patient/'

if len(sys.argv) < 3:
	ofp = open('decipher.txt','w')
else:
	ofp = open(sys.argv[2],'w')
print >> ofp,"Query\tPatient_ID\tLocation\tPhenotypes\tAge_at_Initial_Presentation\tChromosomal_Sex\tMother_is\tFather_is\tPhenotype_count\tVariants_count"
ifp = open(sys.argv[1])
for line in ifp:
	item = line.rstrip().split()
	position = item[0]+':'+item[1]+'-'+item[2]
	values['position']=position

	data = urllib.urlencode(values)
	req = urllib2.Request(ucsc_url, data)
	response = urllib2.urlopen(req)
	with open("ucsc.html", "w") as f:
		    f.write(response.read())

	patient = []
	# open html file
	ifp = open('ucsc.html')
	action = True
	for line in ifp:
		if action and "MAP name='map_side_decipher'" not in line:
			continue
		action = False
		if 'AREA' in line:
			if '/MAP' not in line:
				patient.append(parse_ucsc(line.rstrip()))
			else:
				patient.append(parse_ucsc(line.rstrip()))
				break
	for i,each in enumerate(patient):
		newurl = decipher_url+each[0]
		req = urllib2.Request(newurl)
		response = urllib2.urlopen(req)
		tree = html.fromstring(response.read())
		phenotypes = tree.xpath('//span[@id="phenotype_count"]/text()')
		variants = tree.xpath('//div[@id="variant_count"]/text()')
		overview = tree.xpath('//span[@class="value"]/text()')
		patient[i].extend([item.strip() for item in overview[0:-1]])
		patient[i].extend(phenotypes)
		patient[i].extend(variants)
	for each in patient:
		print >> ofp, position+'\t'+'\t'.join(each)
ifp.close()
ofp.close()
os.remove('ucsc.html')
