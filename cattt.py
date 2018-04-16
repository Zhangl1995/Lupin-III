#coding:utf-8
import numpy as np
from astropy.table import Table, vstack

aor = open('data/ULIRGS_LIRGS_HL.txt','w+')
namelist = open('ID_list/ULIRGS_LIRGS_HL.txt','r')
name = namelist.readlines()
tbl = Table.read('output.tbl', format='ascii')
for each in name:
	dat = open('SHA_data_total/ULIRGS_LIRGS_HL/output'+each.strip()+'.tbl','r+')
	data = dat.readlines()
	if(len(data)>15):
	    dat.close()
	    tbl_i = Table.read('SHA_data_total/ULIRGS_LIRGS_HL/output'+each.strip()+'.tbl', format='ascii')
	    tbl = vstack([tbl,tbl_i], join_type = 'inner')
tbl.write(aor, format='ascii')
namelist.close()
aor.close()
