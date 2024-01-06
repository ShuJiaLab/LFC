# This program gives the expression of the phase modulation of the 
# linear-surface axicons.

# Copyright (c) 2018-2023 by Shawn Hua
# Contact info. xwghua@gmail.com
## --!
import numpy as np
import sys
sys.path.append("../..")
import KymaUti.convMLA.convMLA as convMLA
# import h5py
# import scipy.io as scio

def LensConvex(U0,img_pxsize,wavelen,savefile,datafolder,**param):
	# focal_len):
	# k = 2 * np.pi() / wavelen
	rowsize, colsize = np.shape(U0)
	# widthx = colsize * img_pxsize
	# widthy = rowsize * img_pxsize
	x_range = np.arange(np.ceil(-colsize/2),np.ceil(colsize/2)) * img_pxsize
	y_range = np.arange(np.ceil(-rowsize/2),np.ceil(rowsize/2)) * img_pxsize
	X,Y = np.meshgrid(x_range,y_range)
	# print(x_range,'\n',y_range)

	PhaseOut = np.exp(-1j*np.pi/(wavelen*param['focal_len'])*(X**2 + Y**2))

	if param['MLA']==1:
		print('MLA implemented...')
		import KymaUti.convMLA.convMLA as convMLA
		getPhaseout = convMLA.convMLA(U0 , PhaseOut, img_pxsize, **param)
		PhaseOut = getPhaseout

	Uf = U0 * PhaseOut
	Uyz = Uf[:,int(np.round(colsize/2))]
	print('z -> convex', np.max(np.abs(Uyz)))
	# Uyz = Uyz/np.max(np.abs(Uyz))

	# if savefile:
	# 	scio.savemat(datafolder + '/IND'+'0'*(param['CMPTIND']<10)+str(param['CMPTIND'])+\
	# 		'_convex'+'MLA'*int(param['MLA']) + '.mat', {'Uf':Uf})

	return Uf,Uyz