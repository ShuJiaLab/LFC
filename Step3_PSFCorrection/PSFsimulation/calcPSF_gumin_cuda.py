# import pyopencl as cl
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
#from pyopencl.array import to_device
import numba as nb
import numpy as np
from scipy import integrate
from scipy.special import jn
import time

# import multiprocessing
#from joblib import Parallel,delayed
import os
import sys
sys.path.append("../..")
import KymaUti.kronrod as kr
import KymaUti.gpuConf.getKernel as getKernel

class calcPSF_gumin(object):

	def getPSF(self,p1,p2,p3,fobj,NA,xspace,yspace,wavelen,M,n,n_sample,boundary,calcMode=0):
		# boundary is a limit of centerArea, which represents centerArea(1)
		# print("==Start PSF calculation==")
		self.p1 = p1
		self.p2 = p2
		self.p3 = p3
		self.fobj = fobj
		self.wavelen = wavelen
		self.M = M
		self.boundary = boundary
		self.n = n
		self.n_sample = n_sample
		self.k = 2 * np.pi * n / wavelen
		self.alpha = np.arcsin(min(NA/n,n_sample/n))
		self.xspace = xspace
		self.yspace = yspace
		self.xlength = len(xspace)
		self.ylength = len(yspace)
		self.zeroline = np.zeros((1,int(self.ylength)),dtype=complex)

		self.pattern = np.zeros((self.xlength,self.ylength),dtype=complex)
		self.centerPT = np.ceil(self.xlength/2)
		self.calcMode = calcMode

		self.st = time.perf_counter()


		if calcMode == 0:
			num_cores = multiprocessing.cpu_count()

			self.calcOCTA(boundary)

			self.patternQUAT = self.fmax_complx(self.pattern,np.fliplr(self.pattern))
			self.patternHALF = self.fmax_complx(self.patternQUAT,np.rot90(self.patternQUAT))
			self.PSF_field = self.fmax_complx(self.patternHALF,np.rot90(self.patternHALF,2))

		elif calcMode == 1:
			self.calcOCTA_gpu(boundary)

			self.patternQUAT = self.fmax_complx(self.pattern,np.fliplr(self.pattern))
			self.patternHALF = self.fmax_complx(self.patternQUAT,np.rot90(self.patternQUAT))
			self.PSF_field = self.fmax_complx(self.patternHALF,np.rot90(self.patternHALF,2))

		else:
			print("Calculation Mode error! [ErrInfo] Value our of range.")
			print("[Note] '0' indicates cpu while '1' opencl based gpu.")

	def fmax_complx(self,ma,mb):
		return ((ma+mb)+np.sign(np.abs(ma)-np.abs(mb))*(ma-mb))/2
		
## ===================================================
## sub-functions separated for CPU
	# @nb.jit
	def calcOCTA(self,boundary):
		for a in range(int(boundary)+1,int(self.centerPT+1)):
			x1 = self.xspace[a-1]
			patternLine = self.zeroline
			for b in range(a,int(self.centerPT+1)):
				x2 = self.yspace[b-1]
				#print(x1*1e6,x2*1e6)
				xL2normsq = (((x1+self.M*self.p1)**2+(x2+self.M*self.p2)**2)**0.5)/self.M
				v = self.k*xL2normsq*np.sin(self.alpha)
				u = 4*self.k*(self.p3*1)*(np.sin(self.alpha/2)**2)
				self.u = u
				self.v = v
				Koi = self.M/((self.fobj*self.wavelen)**2)*np.exp(-1j*u/(4*(np.sin(self.alpha/2)**2)))
				U0_re = integrate.quad(self.intgrand_re,0,self.alpha)[0]
				U0_im = integrate.quad(self.intgrand_im,0,self.alpha)[0]
				U0 = U0_re + 1j*U0_im
				#print('b',type(b),'Koi',Koi,'U0',Koi*U0)
				patternLine[0,int(b-1)]=Koi*U0

				N0 = (((len(range(int(boundary)+1,int(self.centerPT+1))))**2)/2)
				N = N0 - len(range(a,int(self.centerPT+1)))**2 /2 + b-a +1
				p = round(N* 100 / N0)
				duration = round(time.perf_counter() - self.st, 2)
				remaining = round(duration * 100 / (0.01 + p) - duration, 2)
				print("Proc:{0}%, Time used:{1}s, Time left:{2}s".format(p, duration, remaining), end="\r")

			self.pattern[int(a-1),:] = patternLine
		print("\n")

	def intgrand_re(self,theta):
		expr_re = (np.sqrt(np.cos(theta)))*(1+np.cos(theta))*\
		(np.cos((self.u/2)*(np.sin(theta/2)**2)/(np.sin(self.alpha/2)**2)))*\
		(jn(0, np.sin(theta)/np.sin(self.alpha)*self.v))*(np.sin(theta))
		return expr_re

	def intgrand_im(self,theta):
		expr_im = (np.sqrt(np.cos(theta)))*(1+np.cos(theta))*\
		(np.sin((self.u/2)*(np.sin(theta/2)**2)/(np.sin(self.alpha/2)**2)))*\
		(jn(0, np.sin(theta)/np.sin(self.alpha)*self.v))*(np.sin(theta))
		return expr_im

## ===================================================
## sub-functions separated for GPU
	def calcOCTA_gpu(self,boundary):
		self.gkx, self.gkw1 = self.calcABSC()

		self.quart_size = int(self.centerPT-boundary)
		self.pattern_quart = np.zeros((self.quart_size*self.quart_size),dtype=np.float64)
		# print('shape of pattern_quart:',np.shape(self.pattern_quart))

		print('Sending data for GPU calculation... <Nvdia GPU Powered> \n')
		if self.n == self.n_sample:
			try:
				pattern_quart_complex = self.gpuCalc('./calcPSFgpu.cu')
			except:
				pattern_quart_complex = self.gpuCalc('./KymaCmpt/PSF/calcPSFgpu.cu')
		else:
			print('*** Applying refractive index mismatch ***')
			try:
				pattern_quart_complex = self.OPTgpuCalc('./KymaCmpt/PSF/OptcalcPSFgpu.cu')
			except:
				pattern_quart_complex = self.OPTgpuCalc('./KymaCmpt/PSF/OptcalcPSFgpu.cu')
		# pattern_quart_complex = self.gpuCalc('./KymaCmpt/PSF/calcPSFgpu.cu')

		self.pattern = np.vstack((self.pattern[0:int(boundary),:],\
			np.hstack((self.pattern[int(boundary):int(self.centerPT),0:int(boundary)],\
				(pattern_quart_complex),\
				self.pattern[int(boundary):int(self.centerPT),int(self.centerPT):])),\
			self.pattern[int(self.centerPT):,:] ))
		
		
		u = 4*self.k*(self.p3*1)*(np.sin(self.alpha/2)**2)/self.n*self.n_sample
		Koi = self.M/((self.fobj*self.wavelen*1e7)**2) *np.exp(-1j*u/(4*(np.sin(self.alpha/2)**2)))
		self.pattern = Koi * self.pattern
		
		# print('========== Execute OpenCL source codes successfully!==========')
		duration = round(time.perf_counter() - self.st, 2)
		print("GPU Powered | Time used: ",duration,"s.\n")


	def calcABSC(self):
		n = 100
		_gkx, _gkw1, _gkw2 = kr.kronrod ( n, 1e-12 )
		
		_gkx = np.hstack((-_gkx[:-1],np.flipud(_gkx)))
		_gkw1 = np.hstack((_gkw1[:-1],np.flipud(_gkw1)))
		#_gkw2 = np.hstack((_gkw2[:-1],np.flipud(_gkw2)))
		#print("_gkx, _gkw1, _gkw2 = ",_gkx)
		#gkx_1, gkw1_1, gkw2_1 = kr.kronrod_adjust(self.alpha/5*0, self.alpha/5*1, 2*n, _gkx, _gkw1, _gkw2)
		#gkx_2, gkw1_2, gkw2_2 = kr.kronrod_adjust(self.alpha/5*1, self.alpha/5*2, 2*n, _gkx, _gkw1, _gkw2)
		#gkx_3, gkw1_3, gkw2_3 = kr.kronrod_adjust(self.alpha/5*2, self.alpha/5*3, 2*n, _gkx, _gkw1, _gkw2)
		#gkx_4, gkw1_4, gkw2_4 = kr.kronrod_adjust(self.alpha/5*3, self.alpha/5*4, 2*n, _gkx, _gkw1, _gkw2)
		#gkx_5, gkw1_5, gkw2_5 = kr.kronrod_adjust(self.alpha/5*4, self.alpha/5*5, 2*n, _gkx, _gkw1, _gkw2)

		#gkx, gkw1, gkw2 = kr.kronrod_adjust(0, self.alpha, 2*n, _gkx, _gkw1, _gkw2)
		gkx = _gkx*self.alpha/2 + self.alpha/2;
		gkw1 = _gkw1*self.alpha/2
		#print(gkx)
		#gkx = np.hstack((gkx_1[:-1],gkx_2[:-1],gkx_3[:-1],gkx_4[:-1],gkx_5))
		#gkw1 = np.hstack((gkw1_1[:-1],gkw1_2[:-1],gkw1_3[:-1],gkw1_4[:-1],gkw1_5))
		return gkx, gkw1


	# def gpuCalc(self,ctx,queue,kernels,boundary):#,
		#dev_p1,dev_p2,dev_p3,
	def gpuCalc(self,cufilename):
		#dev_fobj,dev_k,dev_alpha,dev_M,dev_wavelen,dev_boundary,dev_centerPT,
		#dev_xspace,dev_yspace,dev_pattern_quart_re,dev_pattern_quart_im):
		
		# rowsize,colsize = np.shape(self.pattern_quart)
		matsize = self.quart_size
		blocksize = 16
		gridsize = int(np.ceil(matsize/blocksize))
		
		pattern_quart_re = np.zeros_like(self.pattern_quart).astype(np.float64)
		pattern_quart_im = np.zeros_like(self.pattern_quart).astype(np.float64)

		f = open(cufilename,'r',encoding='utf-8')
		_kernel = ''.join(f.readlines())
		f.close()

		print('*** blocksize, gridsize = ',blocksize,gridsize)

		model = SourceModule(_kernel)
		# print("Kernel compiled successfully...")
		func = model.get_function("calcPSFgpu")
		func(
			cuda.In(np.float64(self.p1*1e6)),cuda.In(np.float64(self.p2*1e6)),cuda.In(np.float64(self.p3*1e6)),
			cuda.In(np.float64(self.fobj)),cuda.In(np.float64(self.k)),cuda.In(np.float64(self.alpha)),
			cuda.In(np.float64(self.M)),cuda.In(np.float64(self.wavelen*1e9)),
			cuda.In(np.int32(self.boundary)),cuda.In(np.int32(self.centerPT)),
			cuda.In((self.xspace*1e6).astype(np.float64)),cuda.In((self.yspace*1e6).astype(np.float64)),
			cuda.In(self.gkx.astype(np.float64)),cuda.In(self.gkw1.astype(np.float64)),
			cuda.InOut(pattern_quart_re),cuda.InOut(pattern_quart_im),
			block=(blocksize,blocksize,1),grid=(gridsize,gridsize,1))


		# prg.calcPSFgpu(gpuKernel.queue,self.pattern_quart.shape,None,
		# 	np.float64(self.p1*1e6),np.float64(self.p2*1e6),np.float64(self.p3*1e6),
		# 	np.float64(self.fobj),np.float64(self.k),np.float64(self.alpha),np.float64(self.M),
		# 	np.float64(self.wavelen*1e9),np.int32(self.boundary),np.int32(self.centerPT),
		# 	gpuKernel.dev_xspace,gpuKernel.dev_yspace,gpuKernel.dev_gkx,gpuKernel.dev_gkw1,
		# 	gpuKernel.dev_pattern_quart_re,gpuKernel.dev_pattern_quart_im)
		# print('Waiting for kernel response ...')
		#evt.wait()
		# pattern_quart_re = np.empty_like(self.pattern_quart)
		# pattern_quart_im = np.empty_like(self.pattern_quart)
		# cl.enqueue_copy(gpuKernel.queue, pattern_quart_re, gpuKernel.dev_pattern_quart_re)
		# cl.enqueue_copy(gpuKernel.queue, pattern_quart_im, gpuKernel.dev_pattern_quart_im)
		pattern_quart_complex = pattern_quart_re.reshape((self.quart_size,self.quart_size)) + \
		1j * pattern_quart_im.reshape((self.quart_size,self.quart_size))

		# pattern_quart_complex = pattern_quart_re + 1j * pattern_quart_im
		#print(pattern_quart_complex,end='\n')
		#print("size",np.shape(pattern_quart_complex))
		#self.pattern_quart.shape)
		return pattern_quart_complex

	def OPTgpuCalc(self,cufilename):
		#dev_fobj,dev_k,dev_alpha,dev_M,dev_wavelen,dev_boundary,dev_centerPT,
		#dev_xspace,dev_yspace,dev_pattern_quart_re,dev_pattern_quart_im):
		
		# rowsize,colsize = np.shape(self.pattern_quart)
		matsize = self.quart_size
		blocksize = 16
		gridsize = int(np.ceil(matsize/blocksize))
		
		pattern_quart_re = np.zeros_like(self.pattern_quart).astype(np.float64)
		pattern_quart_im = np.zeros_like(self.pattern_quart).astype(np.float64)

		f = open(cufilename,'r',encoding='utf-8')
		_kernel = ''.join(f.readlines())
		f.close()

		print('*** blocksize, gridsize = ',blocksize,gridsize)

		model = SourceModule(_kernel)
		# print("Kernel compiled successfully...")
		func = model.get_function("calcPSFgpu")
		func(
			cuda.In(np.float64(self.p1*1e6)),cuda.In(np.float64(self.p2*1e6)),cuda.In(np.float64(self.p3*1e6)),
			cuda.In(np.float64(self.fobj)),cuda.In(np.float64(self.k)),cuda.In(np.float64(self.alpha)),
			cuda.In(np.float64(self.M)),cuda.In(np.float64(self.wavelen*1e9)),
			cuda.In(np.int32(self.boundary)),cuda.In(np.int32(self.centerPT)),
			cuda.In(np.float64(self.n)),cuda.In(np.float64(self.n_sample)),
			cuda.In((self.xspace*1e6).astype(np.float64)),cuda.In((self.yspace*1e6).astype(np.float64)),
			cuda.In(self.gkx.astype(np.float64)),cuda.In(self.gkw1.astype(np.float64)),
			cuda.InOut(pattern_quart_re),cuda.InOut(pattern_quart_im),
			block=(blocksize,blocksize,1),grid=(gridsize,gridsize,1))


		# prg.calcPSFgpu(gpuKernel.queue,self.pattern_quart.shape,None,
		# 	np.float64(self.p1*1e6),np.float64(self.p2*1e6),np.float64(self.p3*1e6),
		# 	np.float64(self.fobj),np.float64(self.k),np.float64(self.alpha),np.float64(self.M),
		# 	np.float64(self.wavelen*1e9),np.int32(self.boundary),np.int32(self.centerPT),
		# 	gpuKernel.dev_xspace,gpuKernel.dev_yspace,gpuKernel.dev_gkx,gpuKernel.dev_gkw1,
		# 	gpuKernel.dev_pattern_quart_re,gpuKernel.dev_pattern_quart_im)
		# print('Waiting for kernel response ...')
		#evt.wait()
		# pattern_quart_re = np.empty_like(self.pattern_quart)
		# pattern_quart_im = np.empty_like(self.pattern_quart)
		# cl.enqueue_copy(gpuKernel.queue, pattern_quart_re, gpuKernel.dev_pattern_quart_re)
		# cl.enqueue_copy(gpuKernel.queue, pattern_quart_im, gpuKernel.dev_pattern_quart_im)
		pattern_quart_complex = pattern_quart_re.reshape((self.quart_size,self.quart_size)) + \
		1j * pattern_quart_im.reshape((self.quart_size,self.quart_size))

		# pattern_quart_complex = pattern_quart_re + 1j * pattern_quart_im
		#print(pattern_quart_complex,end='\n')
		#print("size",np.shape(pattern_quart_complex))
		#self.pattern_quart.shape)
		return pattern_quart_complex


