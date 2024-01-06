from turtle import color
import pyvista as pv
import numpy as np
import imageio as iio
from skimage.transform import resize as imresize
from scipy.ndimage import rotate as imrotate

bloodname = 'red_1000nm_clbu'
zlayerbegin = 35
zlayerend = 55
vol = iio.volread('./beads_selected_1/'+bloodname+'.tif')
vol = np.single(vol[zlayerbegin:zlayerend+1,:,:])
shape0 = vol.shape
vol = imresize(vol,(round(shape0[0]*100/65),shape0[1],shape0[2]))
vol = np.transpose(vol)
rescale_min = -200
rescale_max = 300
vol = (vol-vol.min())/(vol.max()-vol.min())*(rescale_max-rescale_min) + rescale_min
# vol = imrotate(vol,-10,axes=(1,2),reshape=False)
# vol = imrotate(vol,10,axes=(0,2),reshape=False)
# vol = imrotate(vol,140,axes=(0,1),reshape=False)
print(vol.shape)

pv.global_theme.background = 'white'
p = pv.Plotter()
# p.add_mesh(pv.wrap(vol).outline(), color="k",line_width=2)
opacity = [0, 0, 0, 0.1, 0.3, 0.6, 1]
p.add_volume(vol,
    cmap="Reds",
    opacity=opacity,
    show_scalar_bar = False)
imscale_factor = 1.0

p.camera_position = [\
    (-100*imscale_factor, 500.0*imscale_factor, 200.0*imscale_factor), \
    (100, 100, 77),\
    (-0.1332, 0.2632, 0.9427)]
print(p.camera_position)
p.add_bounding_box(line_width=2, color=[0.5,0.5,0.5])
p.show_grid(grid='front',
    location='outer',
    color=[0.85,0.85,0.85],
    # all_edges = True,
    xlabel = '',
    ylabel = '',
    zlabel = '',
    font_size = 0,bold = False,
    axes_ranges = [0,200*0.065,0,200*0.065,0,155*0.065])
p.scale = [1.0,1.0,1.0]
# p.set_position([-100*1.35, 500.0*1.35, 200.0*1.35])
p.window_size = [600,550]
# p.show()
p.show(screenshot='./'+bloodname+'_3d.png')