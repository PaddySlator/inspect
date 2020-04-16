# InSpect toolbox for qMRI analysis #

There are three main commands:
* inspect_seg.m
  * InSpect segmentation algorithm. 
    * Reference: Slator, P. J. et al. InSpect: INtegrated SPECTral Component Estimation and Mapping for Multi-contrast Microstructural MRI. in IPMI 2019 755–766 (Springer, Cham, 2019). https://doi.org/10.1007/978-3-030-20351-1_59 
* inspect_map.m
  * InSpect continuous mapping algorithm 
    * Reference: paper in preparation
* inspect_vox.m
  * Implementation of standard voxelwise spectral estimation and calculation of apparent volumen fractions with spectral integration (i.e. not my algorithm!).
    * Reference (the earliest I can find): English, A.E. et al. Quantitative Two-Dimensional time Correlation Relaxometry.  MRM 22(2), 425–434 (dec 1991). https://doi.org/10.1002/mrm.1910220250

The required inputs are: 
* an image (nifti or matlab array)
* text file of MR acquisition parameters with columns [gradx grady gradz b-value echotime inversiontime repititiontime] (known as "gradechoinv" file)
* binary mask (nifti or matlab array)
* string specifying the kernel to use, e.g. 'DT2' is exp(-bD) * exp(-TE/T2)

For example the following command will run the InSpect segmenetation version on qMRI data "image.nii", which was acquired with the MR parameters in "gradechoinv.txt", on the voxels specified by "mask.nii", using a diffusion-T2 kernel.

`output = inspect_seg(‘image.nii’, ‘gradechoinv.txt’, ‘mask.nii’, ‘DT2’)`

I recommend trying out the `simulated_example_inspect_map` and `simulated_example_inspect_seg` commands to get an idea of how InSpect works. Each command synthesises example datasets then compares the corresponding InSpect version to voxelwise spectral mapping.


Paddy Slator, p.slator@ucl.ac.uk
