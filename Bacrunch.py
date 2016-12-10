#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Biolog.py

Created by Povilas Norvaisas on 2015-02-26.
Copyright (c) 2015. All rights reserved.

"""
try:
	import pip
except ImportError, e:
	print "Module pip not found!"
	print "Please install pip manually to proceed!"
	sys.exit(1)
def install(package):
	pip.main(['install', package])

for mod in ['pip','string','math','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)
		#pass # module doesn't exist, deal with it.
import unicodedata
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc




import scipy.signal as sig
from scipy.fftpack import rfft, irfft, fftfreq
from scipy import interpolate as ip
from collections import OrderedDict
from collections import defaultdict
from collections import Counter
from scipy.optimize import curve_fit
import numpy as np
from operator import itemgetter
from itertools import groupby
import textwrap as tw
import itertools as IT
from multiprocessing import Pool


def numerize(s):
	try:
		if s=='NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s)==0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s


def filename(ifile):
	if ifile.split('.')[0]=='':
		ipat=''
		iname=''
		itype=ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep="\\"
		elif "/" in ifile.split('.')[0]:
			sep="/"
		else:
			ipat=''
			iname=ifile.split('.')[0]
			itype=ifile.split('.')[1]
			return ipat, iname, itype
		allpath=ifile.split('.')[0]
		iname=allpath.split(sep)[-1]
		ipath=allpath.split(sep)[:-1]
		ipat='/'.join(ipath)
		itype=ifile.split('.')[1]
	return ipat, iname, itype

def readxls(ifile,tabnum=0):
	book=xlrd.open_workbook(ifile,formatting_info=False)
	data=book.sheet_by_name([nm for nm in book.sheet_names()][tabnum]) # if 'Magellan' in nm
	sheet=[]
	for r in range(0,data.nrows):
		#if len(data.row_values(r))>200:
			#print set(data.row_values(r))
		sheet.append(data.row_values(r))
	return sheet



def writecsv(data,ofile,delim='\t'):
	f=open(ofile,'wb')
	ofile=csv.writer(f, delimiter=delim) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in data:
		#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()





os.chdir('/Users/Povilas/Projects/2015-Metformin/Worms/4th_bacterial_screen_fill')

ifiles=glob.glob('*.xlsx')
ifiles=[ifl for ifl in ifiles if 'uM' in ifl]

output=[]
ohead=['Plate','Well','Replicate','Drug','OD']
#'Time'
output.append(ohead)

for ifile in ifiles:
	ipat, iname, itype = filename(ifile)
	inms=iname.split()
	plate,rep=inms[3].split('-')
	drug=inms[5]
	#time=inms[3]
	sheet=readxls(ifile)
	print plate,rep,drug#,time
	data=sheet[5:13]
	data=[ [numerize(val) for val in ln if val!=''] for ln in data]# if val!=''
	#print plate,rep,drug

	for r,ln in IT.izip(['A','B','C','D','E','F','G','H'],data):
		for c,val in enumerate(ln):
			well='{}{}'.format(r,c+1)
			output.append([plate,well,rep,drug,val])



writecsv(output,'4th_screen_bacterial_linear_fill.csv',delim=',')



lsy=[['a','1','Cake'],['b','2','Cake'],['a','3','Cake'],['a','1','Cake'],['b','2','Cake']]



ifile='Biolog worm data summary.xlsx'

os.chdir('/Users/Povilas/Projects/2015-Metformin/Worms/Biolog/Worm_data')
ipat, iname, itype = filename(ifile)
sheet=readxls(ifile,tabnum=1)

datacollect=[['Plate','Well','Replicate','Value']]
for row in sheet:
	if 'Rep' in row[0]:
		repl=int(row[0][3:4])
		continue
	elif row[0] in ['PM1','PM2A','PM5']:
		plate=row[0]
		continue
	elif all([val!='' for val in row]):
		for cell,col in IT.izip(row[1:13],range(1,13)):
			well=row[0]+str(col)
			newline=[plate,well,repl,cell]
			datacollect.append(newline)
			print newline


writecsv(datacollect,'Biolog_Worm_linearised_control.csv',delim=',')


os.chdir('/Users/Povilas/Projects/2015-Metformin/Worms/Keio_NGM_full_screen/Rep1')

ifiles=glob.glob('*.xlsx')

output=[]
ohead=['Plate','Well','Replicate','OD']
#'Time'
output.append(ohead)

for ifile in ifiles:
	ipat, iname, itype = filename(ifile)
	inms=iname.split()
	plate=inms[3].replace('Plate','')
	rep=1
	#time=inms[3]
	sheet=readxls(ifile)
	print plate,rep#,time
	data=sheet[5:13]
	#print len(data)
	print 'A1:{}, A12:{}'.format(data[0][0], data[0][11])
	print 'H1:{}, H12:{}'.format(data[7][0], data[7][11])
	data=[ [numerize(val) for val in ln if val!=''] for ln in data]# if val!=''
	#print plate,rep,drug

	for r,ln in IT.izip(['A','B','C','D','E','F','G','H'],data):
		for c,val in enumerate(ln):
			well='{}{}'.format(r,c+1)
			output.append([plate,well,rep,val])



writecsv(output,'NGM_Keio_Rep1.csv',delim=',')




