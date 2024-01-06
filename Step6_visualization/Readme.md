# How to perform 3D visualization

## 1. Software package list

* ImageJ package: ClearVolume 1.4.2
  * For those who have `'const float4 color'` OpenCL bugs in rendering volumes after installing the ClearVolume plugin via FijiImageJ update, here I included an updated jar, so you only need to replace the original .jar file after the installation in the path `<Fiji directory>/jars/clearvolume-1.4.2.jar`.
    Please don't forget to cite the original Github repository and paper: [https://github.com/ClearVolume/clearvolume](https://github.com/xwghua/clearvolume-folked-updated/blob/Patches/url) and **[ClearVolume – Open-source live 3D visualization for light sheet microscopy.](http://www.nature.com/nmeth/journal/v12/n6/full/nmeth.3372.html)** *Loic A. Royer, Martin Weigert, Ulrik Günther, Nicola Maghelli, Florian Jug, Ivo F. Sbalzarini, Eugene W. Myers* , Nature Methods 12, 480–481 (2015) doi:10.1038/nmeth.3372, when you use the software.
* Python package: PyVista
  * Citation: Sullivan et al., (2019). PyVista: 3D plotting and mesh analysis through a streamlined interface for the Visualization Toolkit (VTK). Journal of Open Source Software, 4(37), 1450, https://doi.org/10.21105/joss.01450
* MATLAB 