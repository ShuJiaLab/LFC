import numpy as np
import os
#import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.misc import imresize

from mayavi import mlab
import  moviepy.editor as mpy

import PostOBJ

class makeGIF(object):

	def __init__(self,runtime):
		self.duration = 10
		self.fps = 10

		self.fig = mlab.figure(size=(1700*runtime.scale,1700*runtime.scale),\
			bgcolor=(1,1,1))
		self.PSF_3D_t = lambda t: runtime.PSF_3D[:,:,int(t)]

		self.animation = mpy.VideoClip(self.make_frame, \
			duration = self.duration).resize(1)
		self.animation.write_gif("PSF_3Dmov.gif", fps=self.fps)

	def make_frame(self,t):
		mlab.clf()	
		self.fig.scene._lift()	
		print(round(t,2),np.shape(self.PSF_3D_t(int(t*self.fps)%51)))	
		mlab.surf(self.PSF_3D_t(int(t*self.fps)%51), warp_scale='auto',\
			figure=self.fig)
		return mlab.screenshot(antialiased=True)

class makeSURF(object):

	def __init__(self,runtime,planeInd):
		fig = plt.figure()
		ax = fig.gca(projection='3d')

		surf = ax.plot_surface(runtime.X, runtime.Y, runtime.PSF_3D[:,:,planeInd], \
			cmap=cm.coolwarm,linewidth=0, antialiased=False)

		# ===Customize the z axis.===
		ax.set_zlim(0, 1.01)
		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		# ===Add a color bar which maps values to colors.===
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()

class makeSURF_maya(object):

	def __init__(self, runtime,planeInd):
		s = mlab.mesh(runtime.X, runtime.Y, runtime.PSF_3D[:,:,planeInd]*850*runtime.scale)
		# s = mlab.surf(runtime.PSF_3D[:,:,planeInd], warp_scale='auto')
		mlab.show()


class makeFIO(object):
	"""docstring for makeFIO"""
	def __init__(self, runtime):
		import scipy.io as sio
		sio.savemat('PSF_3D.mat', {'PSF_3D':runtime.PSF_3D})
		
class make3View(object):
	"""docstring for make3View"""
	def __init__(self, runtime,scale,frames,planeInd):
		fig = plt.figure()
		ImgXY = runtime.PSF_3D[:,:,planeInd]
		ImgYZ = runtime.PSF_3D[int(1700*scale/2),:,:]
		ImgYZ = imresize(ImgYZ,(int(1700*scale/2),int(1700*scale)))
		ImgXZ = runtime.PSF_3D[:,int(1700*scale/2),:]
		ImgXZ = imresize(ImgXZ,(int(1700*scale/2),int(1700*scale)))
		plt.subplot(311),plt.imshow(ImgXY)
		plt.subplot(312),plt.imshow(ImgXZ)
		plt.subplot(313),plt.imshow(ImgYZ)
		plt.show()


if __name__ == '__main__':
	print("Beta Model created on 20181012...")
	frames = 200 # better be an even number
	DOF = 4e-6
	scale=0.2
	p3range = [(i-int(frames/2))/int(frames/2)*DOF for i in range(frames+1)]
	runtime = PostOBJ.PostOBJ(scale=scale, p3range=p3range)

	## make a dynamic PSF in a .gif file
	# makeGIF(runtime)

	## make surface for one plane
	# makeSURF(runtime,int(frames/2))

	## make surface/mesh for one plane with mayavi
	# makeSURF_maya(runtime,int(frames/2))

	## make a file output as a .mat format
	# makeFIO(runtime)

	## make three views of the cross section of the PSF
	make3View(runtime,scale,frames,int(frames/2))