#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')

try:
	import pip
except ImportError, e:
	print "Module pip not found!"
	print "Please install pip manually to proceed!"
	sys.exit(1)
def install(package):
	pip.main(['install', package])

for mod in ['pip','string','math','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types','time']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)
		#pass # module doesn't exist, deal with it.


def writecsv(data,ofile,delim='\t'):
	f=open(ofile,'wb')
	ofile=csv.writer(f, delimiter=delim) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in data:
		#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()





import pythoncyc as pc
ecoli=pc.select_organism('ECOLI')

allcofs=ecoli.all_cofactors()

allcofs=[cf for cf in allcofs if 'CPD' in cf]



cfnames=[ecoli[cf].common_name for cf in allcofs]


cofs=['|PYRIDOXAL_PHOSPHATE|']

rxns=ecoli['|PYRIDOXAL_PHOSPHATE|'].cofactors_of
#cmplx=list(set([ecoli[e].enzyme for e in enzs]))

genes=[]

for r in rxns:
	enzyme=ecoli[r].enzyme
	if ecoli[enzyme].components:
		cmps=ecoli[enzyme].components
		#print erxn
		for c in cmps:
			gn=ecoli[c].gene
			print gn
			if isinstance(gn,list):
				for nms in gn:
					gname=ecoli[nms].common_name
					print 'Reaction: {}, Enzyme: {}, Component: {}, Gene: {}'.format(r,enzyme,c,gname)
					genes.append(gname)


cgenes=list(set(genes))

cgenes=[[cg] for cg in cgenes]

writecsv(cgenes,'/Users/Povilas/Projects/2015-Metformin/Worms/Data_final/Genes_using_PLP.csv',delim=',')