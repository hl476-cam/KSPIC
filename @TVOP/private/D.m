function res = D(image,minus1D)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005
% 
% [sx,sy] = size(image);
% 
% Dx = image([2:end,end],:) - image;
% Dy = image(:,[2:end,end]) - image;
% 
% %res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
% res = cat(3,Dx,Dy);

%sssr:
imsz=size(image);

ndim=length(imsz);

if(minus1D)
    switch ndim
        case 1
            res = image([2:end,end]) - image;
            %res=DGradient(image,1);
        case 2
            Dx = image([2:end,end],:) - image;
            %Dx = DGradient(image,1);
%             Dy = image(:,[2:end,end]) - image;
            %Dy = DGradient(image,2);
            res = Dx;%cat(3,Dx,Dy);
        case 3
            Dx = image([2:end,end],:,:) - image;
            %Dx = DGradient(image,1);
            Dy = image(:,[2:end,end],:) - image;
            %Dy = DGradient(image,2);
%              Dz = image(:,:,[2:end,end]) - image;
            %Dz = DGradient(image,3);

            res = cat(4,Dx,Dy);
    end
else
    
    switch ndim
        case 1
            res = image([2:end,end]) - image;
            %res=DGradient(image,1);
        case 2
            Dx = image([2:end,end],:) - image;
            %Dx = DGradient(image,1);
            Dy = image(:,[2:end,end]) - image;
            %Dy = DGradient(image,2);
            res = cat(3,Dx,Dy);
        case 3
            Dx = image([2:end,end],:,:) - image;
            %Dx = DGradient(image,1);
            Dy = image(:,[2:end,end],:) - image;
            %Dy = DGradient(image,2);
             Dz = image(:,:,[2:end,end]) - image;
            %Dz = DGradient(image,3);

            res = cat(4,Dx,Dy,Dz);
    end
end