% KSPIC reconstruction for undersampled subtractive MR angiography data accelerated by compressed sensing,
% parallel imaging and partial Fourier sampling.
% Suitable for data which require the subtraction between two datasets with different contrasts, 
% such as Fresh Blood Imaging (FBI, M. Miyazaki, JMRI, 2000), Flow Sensitive Dephasing 
% (AN Priest, MRM 2012, 2014; Z. Fan, MRM 2009), ASL and CE-MRA.

% The algorithm is described in the paper "Hao Li, et al., Highly Accelerated Subtractive Femoral NCE-MRA 
% using Compressed Sensing with k-space Subtraction, Phase and Intensity Correction, Magn Res Med, 2021".

% Conventional k-space subtraction (KS) and magnitude-subtraction (MS) methods are also added for comparison.

% Please see license files (./licenses) for this program and accompanied codes of:
% 1) mrics.m by Tom Goldstein (TomGoldstein1@gmail.com, https://www.cs.umd.edu/~tomg/projects/split_bregman/
% 2) SPIRiT v0.3 and sparseMRI v0.2 by Michael Lustig (under LICENSE_SPIRiT, mlustig@eecs.berkeley.edu, 
%    https://people.eecs.berkeley.edu/~mlustig/Software.html)
% 3) adapt_array_2d.m by Ricardo Otazo (CBI, New York University)

% Example femoral Fresh Blood Imaging (FBI) MRA datasets from a healthy volunteer are provided in ./data.
% Matrix size: 320(x) * 320(y) * 80(z) * 4(compressed channel)
% Acceleration factor: 20x.

% Copyright, Hao Li, University of Cambridge, United Kingdom. All rights reserved.

clc
clear
close all

addpath(genpath('Utils'));

load(fullfile('data','rawD.mat'));  % Dark-blood FBI dataset
load(fullfile('data','rawB.mat')); % Bright-blood FBI dataset

Param=param_default;

[result,Hdr]=KSPIC(rawB,rawD,Param);


function [result,Hdr]=KSPIC(rawB,rawD,Param)

% Input: 
% rawB, rawD: bright- and dark-blood 3D MRA datasets. 
%             Notice: 1D FFT has been performed along the first dimension (x).
% Param: reconstruction parameters

% Output: 
% result: reconstructed subtracted angiograms 
% Hdr: Information of reconstructed images

%% Acquire the parameters of the data

Hdr.fullSize=size(rawB);    % Data size

pattern = abs(squeeze(rawD(1,:,:,1)));
pattern(pattern>eps) = 1;   
Hdr.pfSize = fun_pfsize(pattern);   % Partial Fourier length in y and z
Hdr.pattern = pattern;              % Sampling pattern
totalNum = sum(pattern(:)~=0);  
Hdr.totalNum = totalNum;            % The number of total sampling points

figure();imshow(pattern');title(['Sampling point number:',num2str(totalNum)]);

calibMask=findFullCenter(pattern);  % Find the fully sampled central region for calibration in SPIRiT
Hdr.caliSize=size(calibMask);       % The size of the calibration region 
Hdr.calibMask=calibMask;            % A mask of the calibration region

Hdr.nCoil=Hdr.fullSize(4);          % Channel number

% Hdr.centreSize=Hdr.caliSize;
Hdr.centreSize=[28,10];             % The size of the central region used for phase correction.


%% Intensity correction


interval=round(Hdr.fullSize(1)/32); % Sampling interval (equally-spaced sampling of axial slices along the x-direction)

Hdr.IC=fun_IC(rawD, rawB, Hdr, Param, interval,true);

%% Reconstruction

gcp;  %Starting parallel pool

timeStart=clock;

fsz=Hdr.fullSize;

[result_KSPIC,result_KS,result_MSIC,result_MS,result_Bright,result_Dark]=deal(zeros(fsz(1),fsz(2),fsz(3)));

sliceD=squeeze(rawD(1,:,:,:));  % Load one slice of dark-blood data for pre-calculating parameters for CS recon

Param=fun_CSprep(sliceD,Param); % Pre-calculate some parameters to reduce computation burden.

fprintf('Reconstruction started.\r');


parfor iSlice=1:fsz(1)
    
    sliceB=squeeze(rawB(iSlice,:,:,:));
    sliceD=squeeze(rawD(iSlice,:,:,:));
    
    result_KSPIC(iSlice,:,:)=fun_recon(sliceD,sliceB,Hdr,Param,'KSPIC');
    
    result_KS(iSlice,:,:)=fun_recon(sliceD,sliceB,Hdr,Param,'KS');
    
    [result_MSIC(iSlice,:,:),result_MS(iSlice,:,:),result_Bright(iSlice,:,:),result_Dark(iSlice,:,:)]=fun_recon(sliceD,sliceB,Hdr,Param,'MS');
    
end

fprintf('Done.\r');

timeEnd=clock;
timeTotal=etime(timeEnd,timeStart);
fprintf('Reconstruction time=%dmin%ds\r',floor(timeTotal/60),round(rem(timeTotal,60)));

%%

result.KS=real(flip(result_KS,1));
result.KS(result.KS<0)=0;
MIP.KS=max(result.KS,[],3);
figure;imshow(MIP.KS,[]);title('K-space subtraction')


result.MS=real(flip(result_MS,1));
result.MS(result.MS<0)=0;
MIP.MS=max(result.MS,[],3);
figure;imshow(MIP.MS,[]);title('Magnitude subtraction')

result.KSPIC=real(flip(result_KSPIC,1));
result.KSPIC(result.KSPIC<0)=0;
MIP.KSPIC=max(result.KSPIC,[],3);
figure;imshow(MIP.KSPIC,[]);title('KSPIC')

save('data\result.mat','result');
save('data\hdr.mat','Hdr');
save('data\Parameter.mat','Param');

end
