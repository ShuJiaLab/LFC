#################################################
# This program aims to generate PSF images
#
# AFTER "OBJECTIVE & TUBE LENS"
#
# This program is an alternative for the MATLAB version
# but with the benifit of OpenCL, we are able to utilize GPU for computation
#
# Developed by Xuanwen Hua (xwghua@gmail.com)
#
################################################
import numpy as np
import sys
sys.path.append("../..")

import KymaUti.gpuConf.getPlatform as getPlatform
_platform = getPlatform.getPlatform()
if (_platform=='N'):
	import KymaCmpt.PSF.calcPSF_gumin_cuda as calcPSF
elif (_platform=='A'):
	import KymaCmpt.PSF.calcPSF_gumin as calcPSF
import KymaCmpt.Free.air as air

import os

os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'

# import h5py

### Set up the scale which determines the sampling frequency and intervals
###
class PostOBJ(object):

	def __init__(self,windowHandle,window_completed,window_Total,\
		NOFILE,xspace=np.linspace(-3000,3000,6000)*3.25e-6,yspace=np.linspace(-3000,3000,6000)*3.25e-6,\
		p3range=np.array([0 + 5e-7]),fobj=0.2/100,NA=1.45,wavelen=680e-9,M=100,n=1.515):

		# self.scale = scale
		# self.x = np.linspace(-850*self.scale,850*self.scale,1700*self.scale)
		# self.y = np.linspace(-850*self.scale,850*self.scale,1700*self.scale)
		self.x = xspace
		self.y = yspace
		self.X,self.Y = np.meshgrid(self.x,self.y)

		self.x1 = self.x
		self.x2 = self.y
		self.p1 = 0
		self.p2 = 0
		# self.p3range = np.array([0 + 5e-7]) 
		self.p3range = p3range
		# self.p3range = np.linspace(-5e-6,5e-6,101) + 5e-7
		# this is tricky where GuMin's model canNOT get phase term at z=0
		# so should be compensated in further transmission

		self.fobj = fobj
		self.NA = NA
		self.wavelen = wavelen
		self.M = M
		self.n = n

		self.boundary = 0#np.ceil(len(self.x)/4)

		self.PSF_3D_0 = np.array([])
		for p3 in self.p3range:
			print('Xuanwen | P3: ',str(p3),', end is 2e-6, interval 1e-7')
			self.PSF_stuff = calcPSF.calcPSF_gumin()
			calcPSF.calcPSF_gumin.getPSF(self.PSF_stuff,
				self.p1,self.p2,p3,
				self.fobj,self.NA,
				self.x1,self.x2,
				self.wavelen,
				self.M,self.n,
				self.boundary,
				calcMode=1)
			#print('PSF_field_out',np.max(np.abs(self.PSF_stuff.PSF_field)))
			self.PSF_intensity = np.abs(self.PSF_stuff.PSF_field)**2
			# print('PSF_intensity',np.max(np.abs(self.PSF_intensity)))
			self.PSF_intensity = self.PSF_intensity/np.max(self.PSF_intensity)
			
			# -!- Here needs the connect the PSF in 3D.
			if len(self.PSF_3D_0)==0:
				self.PSF_3D_0 = self.PSF_stuff.PSF_field
			else:
				self.PSF_3D_0 = np.dstack((self.PSF_3D_0,self.PSF_stuff.PSF_field))

		#print(np.shape(self.PSF_3D[:,85,:]))
		#print(self.PSF_3D)
		#print('PSF_3D',np.max(np.abs(self.PSF_3D)))

		InitParams = {'zdistance':0.005,'stepsize':0.005}
		self.PSF_3D = air.FreeAir(windowHandle,window_completed,window_Total,\
			self.PSF_3D_0,xspace[1]-xspace[0],wavelen,False,NOFILE,**InitParams)[0]
		# self.PSF_3D = self.PSF_3D_0

	# @mlab.animate(delay=200)


if __name__ == '__main__':

	from multiprocessing import Queue
	print("Beta Model created on 20181224...")

	for p3range in [0 + 5e-7]:#np.linspace(-5e-6,5e-6,101) + 5e-7:

		runtime = PostOBJ(windowHandle=Queue(), window_completed=1, window_Total=8,NOFILE='null',\
			p3range=[p3range])

		# from mpl_toolkits.mplot3d.axes3d import Axes3D
		import matplotlib.pyplot as plt
		If = np.absolute(runtime.PSF_3D)**2
		If = If[2901:3101,2901:3101]
		If = If/np.max(If)*255

		fig = plt.figure()
		plt.imshow(If)
		# axes3d = Axes3D(fig)

		# axes3d.plot_surface(runtime.X,runtime.Y,np.abs(runtime.PSF_3D)**2)
		plt.show()

		# import scipy.io as sio
		# sio.savemat('./PSF_FLFM100X/PSF_'+str(np.round((p3range-5e-7)*1e9,2))+'.mat', {'PSF':runtime.PSF_3D})
		




