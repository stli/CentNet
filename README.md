# CentNet
Author: Steffen Limmer steffen.limmer@tu-berlin.de

Matlab code to compute volume and centroid of intersection polytopes for sparse recovery of simplex valued signals. The code was used to simulate the proposed Laplace transform techniques. Please cite the paper below whenever you use this code. Note: numerical instabilities may occur for larger input dimensions.

@article{LiSt17, author = {Limmer, S. and Stanczak, S.}, title = {Optimal deep neural networks for sparse recovery via Laplace techniques}, year = {2017}, journal={arXiv preprint}
}

[General Information]

    This software was tested on MATLAB Version R2016a under Win7 64bit
    Integration by simplicial decomposition is build based on the following functions:
    -the functions 'lcon2vert.m', 'qlcon2vert.m' and 'vert2lcon.m' are Copyright (c) 2015, Matt Jacobson All rights reserved (see license_mj.txt).
    -the function 'gamrnd_m' and 'gen_vec.m' are from the Randomized Algorithms Control Toolbox (RACT), all credit goes to the original authors
    -the function qhullmx comes with Matlab and must be copied to the working directory from matlabroot\toolbox\matlab\polyfun\private
    MATLAB (R) is a registered trademark of The MathWorks, Inc.
    Last updated: 31.08.2017

[Disclaimer] The code is provided "as is", without warranty of any kind. In no event shall the authors or copyright holders be liable for any claim, damages or other liability. You can redistribute it (except for 'lcon2vert.m', 'qlcon2vert.m', 'vert2lcon.m', 'gamrnd_m', 'gen_vec.m', qhullmx) and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
