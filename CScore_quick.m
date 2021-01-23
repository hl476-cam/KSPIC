% A quick version of CS reconstruction with 2 iterations and no parallel imaging reconstruction.
% imaging reconstruction is not added.
% Developed based on 
% 1) mrics.m by Tom Goldstein (TomGoldstein1@gmail.com),
% 2) sparseMRI by Michael Lustig (mlustig@eecs.berkeley.edu, https://people.eecs.berkeley.edu/~mlustig/Software.html)
% 3) Combined MS-based CS recon by Stanislas Rapacchi (S. Rapacchi, MRM, 2014)


function [u,v] = CScore_quick(f, g)



%% Parameters

lambda=4e-4; % TV
nu=1e-9;      % Feedback
threshTV=0.8e-1;

%% Preparations
% 2D FFT
kspace_to_imspace2D = @(x) fftshift(fftshift(flip(flip(ifft2(fftshift(fftshift(x,1),2)),1),2),1),2);
imspace_to_kspace2D = @(x) ifftshift(ifftshift(fft2(flip(flip(ifftshift(ifftshift(x,1),2),1),2)),1),2);

fsz=size(f);
gsz=size(g);

% Reserve memory for the auxillary variables
u = zeros(fsz);
v = zeros(gsz);
% A = zeros(fsz);
% B = zeros(fsz);
D = zeros(fsz);

%% Get sampling pattern
R=zeros(fsz);
R(logical(abs(f)>eps))=1;
acqlines = squeeze(sum(R,3));
R([acqlines>=1])=1;

S=zeros(gsz);
S(logical(abs(g)>eps))=1;
acqlines = squeeze(sum(S,3));
S([acqlines>=1])=1;

%% Normalization
rf = kspace_to_imspace2D(R.*f);
rg = kspace_to_imspace2D(S.*g);

normFactor = 1e0/max(abs(rf(:)));
normFactor = min(normFactor,1e0/max(abs(rg(:))));
f = normFactor*f;
f0 = f;
g = normFactor*g;
g0 = g;

% update
rf = normFactor*rf;
rg = normFactor*rg;

%% Build kernel
uker = zeros(fsz);
uker(round(fsz(1)/2),round(fsz(2)/2),:)=4;
uker(round(fsz(1)/2)-1,round(fsz(2)/2),:)=-1;
uker(round(fsz(1)/2)+1,round(fsz(2)/2),:)=-1;
uker(round(fsz(1)/2),round(fsz(2)/2)-1,:)=-1;
uker(round(fsz(1)/2),round(fsz(2)/2)+1,:)=-1;

iuker = imspace_to_kspace2D(uker);

ukerf = R+(lambda)*iuker+nu;
ukerg = S+(lambda)*iuker+nu;

TV= TVOP(1);    % Total variation

x = zeros([fsz,2]);
bx = zeros([fsz,2]);
y = zeros([gsz,2]);
by = zeros([gsz,2]);

%%  Split-Bregman reconstruction
A  = rf + lambda*(TV'.*(x-bx)) + nu*u ;
B = rg + lambda*(TV'.*(y-by))  + nu*v ;


u = kspace_to_imspace2D(imspace_to_kspace2D(A)./ukerf);
v = kspace_to_imspace2D(imspace_to_kspace2D(B)./ukerg);

% update z
dx = TV*u;
x = shrinkTV( dx+bx,threshTV);
dy = TV*v;
y = shrinkTV( dy+by,threshTV);

bx = bx+dx-x;
by = by+dy-y;

f = f+f0-R.*imspace_to_kspace2D(u);
g = g+g0-S.*imspace_to_kspace2D(v);

rf = kspace_to_imspace2D(R.*f);
rg = kspace_to_imspace2D(S.*g);

A  = rf + lambda*(TV'.*(x-bx)) + nu*u ;
B = rg + lambda*(TV'.*(y-by))  + nu*v ;

u = kspace_to_imspace2D(imspace_to_kspace2D(A)./ukerf);
v = kspace_to_imspace2D(imspace_to_kspace2D(B)./ukerg);

%%%%%%%%% Another iteration %%%%%%%%%%%%%%%%%%%
% update z
dx = TV*u;
x = shrinkTV( dx+bx,threshTV);
dy = TV*v;
y = shrinkTV( dy+by,threshTV);

bx = bx+dx-x;
by = by+dy-y;

f = f+f0-R.*imspace_to_kspace2D(u);
g = g+g0-S.*imspace_to_kspace2D(v);

rf = kspace_to_imspace2D(R.*f);
rg = kspace_to_imspace2D(S.*g);

A  = rf + lambda*(TV'.*(x-bx)) + nu*u ;
B = rg + lambda*(TV'.*(y-by))  + nu*v ;

u = kspace_to_imspace2D(imspace_to_kspace2D(A)./ukerf);
v = kspace_to_imspace2D(imspace_to_kspace2D(B)./ukerg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% undo the normalization so that results are scaled properly
u = u/normFactor;
v = v/normFactor;

v=fftshift(v,2);
u=fftshift(u,2);


return;


function [zs] = shrinkTV(z,lambda)

s=0;
for i=1:size(z,4)
    s = s + z(:,:,:,i).*conj(z(:,:,:,i));
end
s = sqrt(s);
ss = s-lambda;
ss = ss.*(abs(s)>lambda);

s = s+(s<lambda);
ss = ss./s;

zs = repmat(ss,[1 1 1 size(z,4)]).*z;

return;
