% CS reconstruction core using Split-Bregman and SPIRiT.
% Developed based on 
% 1) mrics.m by Tom Goldstein (TomGoldstein1@gmail.com),
% 2) SPIRiT and sparseMRI by Michael Lustig (mlustig@eecs.berkeley.edu, https://people.eecs.berkeley.edu/~mlustig/Software.html)
% 3) Combined MS-based CS recon by Stanislas Rapacchi (S. Rapacchi, MRM, 2014)


function u = CScore_normal(f,Param,CSReconPara)

S=Param.S;
kCalibMask=Param.kCalibMask;
uker=CSReconPara.uker;

lambda=CSReconPara.lambda;
nu=CSReconPara.nu;
threshTV=CSReconPara.threshTV;

fsz=size(f);    %Size
Sneg=repmat(ones(size(S(:,:,1)))-S(:,:,1),[1 1 fsz(3)]);

kspace_to_imspace2D = @(x) fftshift(fftshift(flip(flip(ifft2(fftshift(fftshift(x,1),2)),1),2),1),2);
imspace_to_kspace2D = @(x) ifftshift(ifftshift(fft2(flip(flip(ifftshift(ifftshift(x,1),2),1),2)),1),2);

% Reserve memory for the auxillary variables
u = zeros(fsz);

%% Normalization

rf = kspace_to_imspace2D(S.*f);
normFactor = 1e0/max(abs(rf(:)));
f = normFactor*f;
f0 = f;
rf = normFactor*rf;

%% SPIRIT kernel
kSize = [5,5];
Ncoils = fsz(3);
kernel = zeros([kSize,Ncoils,Ncoils]);

kCalib = mycrop(f,size(kCalibMask),'c');
[AtAG,] = dat2AtA(kCalib,kSize);

for n=1:Ncoils
    kernel(:,:,:,n) = calibrate(AtAG,kSize,Ncoils,n,2);
end

GOP = SPIRiT(kernel, 'conv',fsz(1:2));

%% Total Variation
TV= TVOP(1);

y = zeros([fsz,2]);
b = zeros([fsz, 2]);

%%  Split-Bregman reconstruction with POCS SPIRiT

for outer = 1:Param.nBreg
    for inner = 1:Param.nInner
        rhs = rf + lambda*(TV'.*(y-b))+ nu*u ;
        ku= imspace_to_kspace2D(rhs)./uker;
        u = kspace_to_imspace2D(Sneg.*(ku+GOP.*ku)+f0);             
        
        d = TV*u;
        y = shrinkTV( d+b,threshTV);
        b = b+d-y;      
        
        %         figure(99);imshow(abs(squeeze(v(:,:,1))),[]);
    end
    f = f+f0-S.*imspace_to_kspace2D(u);
    rf = kspace_to_imspace2D(S.*f);
end

u = u/normFactor;

% ku=ku/normFactor;

u=fftshift(u,2);

return;



function [as] = shrinkL1aL2(a,thresha,alpha)
% p=1;
l1smooth =0;%1e-12;
p=1;
sl1 = p*abs(a).*(a.*conj(a)+l1smooth).^((p-1)/2);
p=2;
sl2 = p*abs(a).*(a.*conj(a)+l1smooth).^((p-1)/2);

s=sl1+alpha*sl2;

%%figure, hist(s(:),100);

ss = s-thresha;
ss = ss.*(ss>0);

s = s+(s<thresha);
ss = ss./s;

as = ss.*a;

function [zs] = shrinkTV(z,lambda)

s=0;
for i=1:size(z,4)
    s = s + z(:,:,:,i).*conj(z(:,:,:,i));
end
s = sqrt(s);

%%figure, hist(s(:),100);

ss = s-lambda;
ss = ss.*(abs(s)>lambda);

s = s+(s<lambda);
ss = ss./s;

zs = repmat(ss,[1 1 1 size(z,4)]).*z;

return;
