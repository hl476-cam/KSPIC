function [kspace_lp,kspace_hp] = fun_homo_filter(partial_kspace,Hdr,filterLength)
% Homodyne high-pass and low-pass filter for data with partial Fourier in two directions

fullRow=Hdr.fullSize(2);    % Full size of data matrix
fullCol=Hdr.fullSize(3);

partRow=Hdr.pfSize(1); % Partial sampled region size
partCol=Hdr.pfSize(2);

centreRow=Hdr.centreSize(1);    % The size of the central region for phase correction
centreCol=Hdr.centreSize(2);

blankRow=fullRow-partRow;   % Blank region lengths
blankCol=fullCol-partCol;

filterLengthRow=filterLength;   % Lenth of hanning filter
filterLengthCol=filterLength;

%% Low-pass filter
mask_lp_init=zeros(fullRow,fullCol);

mask_lp_init(round((fullRow-centreRow)/2+1):round((fullRow+centreRow)/2),round((fullCol-centreCol)/2+1):round((fullCol+centreCol)/2))=1;

hanFilter=hanning(filterLengthRow)*hanning(filterLengthCol)';   % Apply a 2D hanning filter to smooth the low-pass mask
mask_lp_filtered = imfilter(mask_lp_init,hanFilter,'same');
mask_lp_filtered=mask_lp_filtered/max(mask_lp_filtered(:));

mask_lp=mask_lp_filtered;

%% Hight-pass filter (Nine blocks)
mask_hp_init=zeros(fullRow,fullCol);
mask_hp_init(1:blankRow,1:blankCol)=2; % 1 top-left
mask_hp_init(1:blankRow,blankCol+1:partCol)=2; %2
mask_hp_init(1:blankRow,partCol+1:fullCol)=0; %3

mask_hp_init(blankRow+1:partRow,1:blankCol)=2; %4
mask_hp_init(blankRow+1:partRow,blankCol+1:partCol)=1; %5 centre
mask_hp_init(blankRow+1:partRow,partCol+1:fullCol)=0; %6

mask_hp_init(partRow+1:fullRow,1:blankCol)=0; %7
mask_hp_init(partRow+1:fullRow,blankCol+1:partCol)=0; %8
mask_hp_init(partRow+1:fullRow,partCol+1:fullCol)=0; %9

mask_hp=mask_hp_init;

%% Do filtering
kspace_hp = mask_hp .* partial_kspace;
kspace_lp = mask_lp .* partial_kspace;


%% 

% figure;
% subplot(3,2,1);imshow(abs(partial_kspace').^0.1,[]);title('original k-space');
% subplot(3,2,2);imshow(mask_lp_init',[]);title('Inital low-pass mask');
% subplot(3,2,3);imshow(mask_lp_filtered',[]);title('Filtered low-pass mask');
% subplot(3,2,4);imshow(mask_hp_init',[]);title('High-pass mask');
% subplot(3,2,5);imshow(abs(kspace_hp').^0.1,[]);title('High-pass filtered k-space');
% subplot(3,2,6);imshow(abs(kspace_lp').^0.1,[]);title('Low-pass filtered k-space');


end

