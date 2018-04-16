#coding:utf-8
from astropy.io import ascii
import numpy as np

name = 'AGN_QSO_Rad'
Specfile = ascii.read('data_IRS_stare/'+name+'.txt')
dat = open('data_IRS_stare/aaAOR_num.txt','a+')
num = len(Specfile['AORKEY'])
dat.write('%-20s\n'%name)
dat.write('----------------------------------------\n')
while(len(Specfile['AORKEY']) != 0):
	m = len(Specfile['AORKEY'])
	AORKEY = Specfile['AORKEY'][0]
	mask = (Specfile['AORKEY'] != AORKEY)
	Specfile = Specfile[mask]
	n = m - len(Specfile[['AORKEY']]) - 1
	num = num - n
	dat.write('%-25s  %-5s\n'%(AORKEY,np.str(n+1)))
dat.write('----------------------------------------\n')
dat.write('%-23s  %-5s\n'%('AORKEY_num',np.str(num)))
dat.write('----------------------------------------\n\n\n')
dat.close()
