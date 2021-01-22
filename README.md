# KSPIC
KSPIC (k-space subtraction with intensity and phase correction) reconstruction for undersampled subtractive MR angiography data accelerated by compressed sensing, parallel imaging and partial Fourier sampling.

Suitable for data which require the subtraction between two datasets with different contrasts, such as Fresh Blood Imaging (FBI, M. Miyazaki, JMRI, 2000), Flow Sensitive Dephasing (AN Priest, MRM 2012, 2014), ASL and CE-MRA.

Split-Bregman and POCS-SPIRiT were used for CS and PI reconstruction respectively. Conventional k-space subtraction and magnitude-subtraction methods are also added for comparison.

The algorithm is described in the paper "Highly Accelerated Subtractive Femoral NCE-MRA using Compressed Sensing with k-space Subtraction, Phase and Intensity Correction, Magn Res Med, 2021".
