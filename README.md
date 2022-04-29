# multi-echo-qsm
This repository contains the Matlab code used to run the analyses in:

Emma Biondetti, Anita Karsa, Francesco Grussu, Marco Battiston, Marios C. Yiannakas, David L. Thomas, Karin Shmueli, "Multi-echo Quantitative Susceptibility Mapping: How to Combine Echoes for Accuracy and Precision at 3 T"
bioRxiv 2021.06.14.448385; doi: https://doi.org/10.1101/2021.06.14.448385 

## Dependencies
The following toolboxes are needed (or recommended) to run this code:
* MEDI toolbox (required): http://pre.weill.cornell.edu/mri/pages/qsm.html
* Tools for NIfTI and ANALYZE image (optional, used for reading/writing NIfTI images): https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
* Image Processing Toolbox for MATLAB (optional, used for image resampling)
* Statistics and Machine Learning Toolbox for MATLAB (optional, used for calculating normal distributions)
* Parallel Computing Toolbox for MATLAB (optional, used to accelerate the calculation of R2* maps by nonlinear fitting)
