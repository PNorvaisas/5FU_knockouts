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

def readxls(ifile):
	book=xlrd.open_workbook(ifile,formatting_info=False)
	data=book.sheet_by_name([nm for nm in book.sheet_names() if 'Magellan' in nm][0])
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


