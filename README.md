# InSpect toolbox for qMRI analysis #

Author: Paddy Slator (p.slator@ucl.ac.uk), Centre for Medical Image Computing, University College London 

## Required Matlab toolboxes ##
* Optimisation Toolbox
* Deep Learning Toolbox
* Image Processing Toolbox 

## Main Functions ##
There are three main commands:

* inspect_seg.m
  * InSpect segmentation algorithm. 
    * Reference: Slator, P. J. et al. InSpect: INtegrated SPECTral Component Estimation and Mapping for Multi-contrast Microstructural MRI. in IPMI 2019 755–766 (Springer, Cham, 2019). https://doi.org/10.1007/978-3-030-20351-1_59 
* inspect_map.m
  * InSpect continuous mapping algorithm 
    * Reference: Slator, P. J. et al. Data-Driven Multi-Contrast Spectral Microstructure Imaging with InSpect: INtegrated SPECTral Component Estimation and Mapping. Medical Image Analysis (2021). https://doi.org/10.1016/j.media.2021.102045
* inspect_vox.m
  * Implementation of standard voxelwise spectral estimation and calculation of apparent volume fractions with spectral integration (i.e. not my algorithm!).
    * Reference (the earliest I can find): English, A.E. et al. Quantitative Two-Dimensional time Correlation Relaxometry.  MRM 22(2), 425–434 (dec 1991). https://doi.org/10.1002/mrm.1910220250

The required inputs are: 
* an image (nifti or matlab array)
* text file of MR acquisition parameters with columns [gradx grady gradz b-value echotime inversiontime repititiontime] (known as "gradechoinv" file, see section below for more details)
* binary mask (nifti or matlab array)
* string specifying the kernel to use, e.g. 'DT2' is exp(-bD) * exp(-TE/T2)

For example the following command will run the InSpect segmenetation version on qMRI data "image.nii", which was acquired with the MR parameters in "gradechoinv.txt", on the voxels specified by "mask.nii", using a diffusion-T2 kernel.

`output = inspect_seg(‘image.nii’, ‘gradechoinv.txt’, ‘mask.nii’, ‘DT2’)`

## Grad Echo Inv Files ##
The gradechoinv file encodes the MR acquisition parameters for each imaging volume. The format is essentially an MRtrix-style grad file with extra columns. Each row corresponds to one imaging volume, and contains space-separated values [x y z b TE TI TR], where 
* x, y, z are the diffusion gradient directions
* b is the b-value
* TE is the echo time 
* TI is the inversion time
* TR is the repitition time

If one of these contrast mechanisms isn't varied in an experiment, you can pad the relevent column with NaN's. If the TI isn't varied then the final two columns can be removed, for example for a T2-diffusion acquisition the rows are [x y z b TE]. There are some example gradechoinv files in the "examples" folder.

## Examples ##
The examples folder also contains demonstrations of the functions. E.g. the `simulated_example_inspect_map` and `simulated_example_inspect_seg` commands synthesise example datasets, and then compare the corresponding InSpect version to voxelwise spectral mapping.

<!--## Adding new kernels ##
The following kernels are currently implemented:

* 

To add a new kernel to the repo you'll need to create/update the following files:

* kernels/Kernel["name\_of\_new\_kernel"].m: create a function that takes a gradechoinv file and model parameters as input and outputs the signal for the new kernel.
* utilities/GetKernelParameterStrings.m: add strings for all parameters of the new kernel - making sure to match existing parameters names 
* update the list above in README.md!-->

## Noise estimation ##
There are multiple options for estimating the noise:

* Estimate the noise in each voxel by calculating the standard deviation of the measurements with highest signal, e.g. b=0 for diffusion, minimum TE for T2 relaxometry, and maximum TI for T1 inversion recovery [default option].
* Use a fixed sigma/SNR value by setting `options.SNR='fixed'` and `options.fixedsigma` or `options.fixedSNR` to the desired sigma/SNR value respectively
* Use a voxelwise noise map, e.g. the output from MP-PCA denoising, by setting: `options.fixedsigma = noisemap_path.nii.gz`


## Citations ##
If you use the code then please cite the appropriate paper:

* *inspect_map.m*: Slator, P. J. et al. *Data-Driven Multi-Contrast Spectral Microstructure Imaging with InSpect: INtegrated SPECTral Component Estimation and Mapping.* **Medical Image Analysis** (2021). https://doi.org/10.1016/j.media.2021.102045
* *inspect_seg.m*: Slator, P. J. et al. *InSpect: INtegrated SPECTral Component Estimation and Mapping for Multi-contrast Microstructural MRI.* in **IPMI** 2019 755–766 (Springer, Cham, 2019). https://doi.org/10.1007/978-3-030-20351-1_59 

* *inspect_vox.m*: (not my algorithm - this is the earliest reference I can find!): English, A.E. et al. *Quantitative Two-Dimensional time Correlation Relaxometry.*  **MRM** 22(2), 425–434 (dec 1991). https://doi.org/10.1002/mrm.1910220250

## Acknowledgements ##
Mark Does's MERA toolbox (https://github.com/markdoes/MERA) was very helpful for understanding and testing when developing this code. 

