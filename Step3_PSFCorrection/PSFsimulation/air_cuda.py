# Kyma phase file for Fresnel Transmission

import numpy as np
import cupy as cp

def FreeAir(U0,img_pxsize,wavelen,**param):
	# zdistance,stepsize):
	k = 2.0 * np.pi / wavelen
	rowsize, colsize = np.shape(U0)
	widthx = colsize * img_pxsize
	widthy = rowsize * img_pxsize
	widthx_range = np.arange(np.ceil(-colsize/2),np.ceil(colsize/2))
	widthy_range = np.arange(np.ceil(-rowsize/2),np.ceil(rowsize/2))
	u, v = np.meshgrid(widthx_range,widthy_range)

	Uf = U0
	Uyz_slice = Uf[:,int(np.round(colsize/2))] 
	# Uyz_slice = Uf[:,np.where(np.abs(Uf)==np.max(np.abs(Uf)))[1][0]]

	Uyz_slice = Uyz_slice/np.max(np.abs(Uyz_slice))
	Uyz = Uyz_slice

	Chirp = np.exp(1j*k*param['stepsize']) * \
	np.exp( -(1j*np.pi)*(wavelen*param['stepsize'])*((u/widthx)**2 + (v/widthy)**2) )

	Chirp_gpu = cp.asarray(Chirp)
	Uf_gpu = cp.asarray(Uf)
	Uyz_gpu = cp.asarray(Uyz)
	Uyz_slice_gpu = cp.asarray(Uyz_slice)

	for z in np.arange(np.round(param['zdistance']/param['stepsize'])):
		# Uf = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(\
		# 	np.fft.fftshift(np.fft.fft2(np.fft.fftshift(Uf)))*Chirp)))
		
		Uf_gpu = cp.fft.fftshift(cp.fft.ifft2(cp.fft.fftshift(cp.fft.fftshift(cp.fft.fft2(cp.fft.fftshift(Uf_gpu)))*Chirp_gpu)))

		Uyz_slice_gpu = Uf_gpu[:,int(np.round(colsize/2))]
		# Uyz_slice = Uf[:,np.where(np.abs(Uf)==np.max(np.abs(Uf)))[1][0]]

		print('(z = ',z,')\t\t max value=',cp.max(cp.absolute(Uyz_slice_gpu)),'\t[CUDA powered]')
		# Uyz_slice = Uyz_slice/np.max(np.abs(Uyz_slice))
		# print(np.max(np.abs(Uyz_slice)))
		Uyz_gpu = cp.vstack((Uyz_gpu,Uyz_slice_gpu))
	Uyz_gpu = Uyz_gpu[1:,:]

	return Uf_gpu.get(),Uyz_gpu.get()

