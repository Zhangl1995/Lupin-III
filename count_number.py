#coding:utf-8
from astropy.io import ascii
import numpy as np

name = 'Local_Group_Gal'
Specfile = ascii.read('data_IRS/'+name+'.txt')
dat = open('data_IRS/aacount_num.txt','a+')
num = len(Specfile['AORKEY'])
dat.write('%-20s\n'%name)
dat.write('%-23s  %-5s\n'%('AOR_num',np.str(num)))
dat.write('----------------------------------------\n')
while(len(Specfile['AORKEY']) != 0):
	m = len(Specfile['AORKEY'])
	Target = Specfile['Target_name'][0]
	mask = (Specfile['Target_name'] != Target)
	Specfile = Specfile[mask]
	n = m - len(Specfile[['AORKEY']]) - 1
	num = num - n
	if(n>0):
		dat.write('%-25s  %-5s\n'%(Target,np.str(n+1)))
dat.write('----------------------------------------\n')
dat.write('%-23s  %-5s\n'%('Target_num',np.str(num)))
dat.write('----------------------------------------\n\n\n')
dat.close()
