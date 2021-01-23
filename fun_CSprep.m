function Param=fun_CSprep(f,Param)
% Pre-calculate parameters to slightly accelerate the reconstruction

fsz=size(f);    %Size
S=zeros(fsz);
S(logical(abs(f)>eps))=1;
acqlines = squeeze(sum(S,3));
S(acqlines>=1)=1;   % Sampling mask

imspace_to_kspace2D = @(x) ifftshift(ifftshift(fft2(flip(flip(ifftshift(ifftshift(x,1),2),1),2)),1),2);

uker = zeros(fsz);
uker(round(fsz(1)/2),round(fsz(2)/2),:)=4;
uker(round(fsz(1)/2)-1,round(fsz(2)/2),:)=-1;
uker(round(fsz(1)/2)+1,round(fsz(2)/2),:)=-1;
uker(round(fsz(1)/2),round(fsz(2)/2)-1,:)=-1;
uker(round(fsz(1)/2),round(fsz(2)/2)+1,:)=-1;

kCalibMask=findFullCenter(S(:,:,1));

iuker = imspace_to_kspace2D(uker);

% ukergSub = S+(Param.sub.lambda)*iuker+Param.sub.nu+Param.sub.mu;
ukerSub = S+(Param.sub.lambda)*iuker+Param.sub.nu;
ukerMag = S+(Param.mag.lambda)*iuker+Param.mag.nu;

Param.S=S;
Param.kCalibMask=kCalibMask;
Param.sub.uker=ukerSub;
Param.mag.uker=ukerMag;

end

