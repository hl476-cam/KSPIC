# KSPIC
This is a Matlab implementation of KSPIC (k-space subtraction with intensity and phase correction) reconstruction for undersampled subtractive MR angiography data accelerated by compressed sensing, parallel imaging and partial Fourier sampling.

Suitable for data which require the subtraction between two datasets with different contrasts, such as Fresh Blood Imaging (FBI, M. Miyazaki, JMRI, 2000), Flow Sensitive Dephasing (AN Priest, MRM 2012, 2014), ASL and CE-MRA.

Split-Bregman and POCS-SPIRiT are used for CS and PI reconstruction respectively. 

Conventional k-space subtraction and magnitude-subtraction methods are also added for comparison.

The algorithm is described in the paper "Hao Li, et al., Highly Accelerated Subtractive Femoral NCE-MRA using Compressed Sensing with k-space Subtraction, Phase and Intensity Correction, Magn Res Med, 2021".


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Please see license files (./licenses) for this program and accompanied codes (which also include own modifications) of:
1) mrics.m by Tom Goldstein (TomGoldstein1@gmail.com, https://www.cs.umd.edu/~tomg/projects/split_bregman/
2) SPIRiT v0.3 and sparseMRI v0.2 by Michael Lustig (under LICENSE_SPIRiT, mlustig@eecs.berkeley.edu, 
   https://people.eecs.berkeley.edu/~mlustig/Software.html)
3) adapt_array_2d.m by Ricardo Otazo (CBI, New York University)

Example femoral Fresh Blood Imaging (FBI) MRA datasets from a healthy volunteer are provided in ./data.
Matrix size: 320(x) * 320(y) * 80(z) * 4(compressed channel)
Acceleration factor: 20x.

Copyright, Hao Li, University of Cambridge, United Kingdom. All rights reserved.
