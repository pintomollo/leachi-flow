29 March 2013, author: Filippo Piccinini (f.piccinini@unibo.it).

NOTES BEFORE BUILDING THE MOSAICS:

To build the mosaic:
- Read the file “_license.txt”.
- Open, read and change the Matlab file named: "ParametersByUser.m".
- Open, read and play the Matlab file named: "START.m".
Requirements: MATLAB 2009b and Image Processing Toolbox 8.1 or later versions. MATLAB must be installed at 32-bit.  

The input video must be recorded with no compression and no interlacing.
The final mosaic will be saved in the "OUTPUT" folder.
In order to perform the flat-field correction, a 2D vignetting surface must be saved as a Matlab matrix as "NAME.m" and stored in the "VIGNETTINGFUNCTION" folder. The dimension of the matrix must be the same as the images to be corrected.
In order to perform the White Balancing, a RGB image (3 channels) referring to an empty field (image acquired without placing objects or placing an homogeneous reference object on the sample holder of the microscope) must be copied in the "WHITEBALANCING" folder. 

The functions used to register the images (that is, "STFeatureExtract.mexw32" and "mexFunctionLKT.mexw32") have been built by importing functions of the C++ OpenCV (see file "_licenseOpenCV.txt") and to work they simply require that the files with extension ".dll" are already present in the main folder.

We thank Peter D. Kovesi that with his Matlab functions (MATLAB and Octave Functions for Computer Vision and Image Processing, freely available in the website: http://www.csse.uwa.edu.au/~pk/research/matlabfns/, see "_licenseMIT.txt") has contributed and inspired many functions of this tool.





