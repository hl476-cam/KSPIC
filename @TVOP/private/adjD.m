function res = adjD(y,minus1D)

% res = zeros(size(y,1),size(y,2));
res = zeros(size(y));

%y1 = ones(imsize)*y(1)/sqrt(prod(imsize));
%yx = (reshape(y(2:prod(imsize)+1), imsize(1), imsize(2)));
%yy = (reshape(y(prod(imsize)+2:end), imsize(1), imsize(2)));

% % res = adjDx(y(:,:,1)) + adjDy(y(:,:,2));


%sssr:
imsz=size(y);

ndim=length(imsz);
if(minus1D)
     
    switch ndim
        case [1,2]
            res = y([1,1:end-1]) - y;
            res(1) = -y(1);
            res(end) = y(end-1);
        case 3
            disp('TV: adjD: should never happens');
%             res = adjD2Dx(y(:,:,1)) + adjD2Dy(y(:,:,2));
        case 4
          res = adjD3Dx(y(:,:,:,1)) + adjD3Dy(y(:,:,:,2));    
    end
else
    
    switch ndim
        case [1,2]
            res = y([1,1:end-1]) - y;
            res(1) = -y(1);
            res(end) = y(end-1);
        case 3
            res = adjD2Dx(y(:,:,1)) + adjD2Dy(y(:,:,2));
        case 4
          res = adjD3Dx(y(:,:,:,1)) + adjD3Dy(y(:,:,:,2))+ adjD3Dz(y(:,:,:,3));    
    end
end

return;


function res = adjD2Dy(x)
res = x(:,[1,1:end-1]) - x;
res(:,1) = -x(:,1);
res(:,end) = x(:,end-1);

function res = adjD2Dx(x)
res = x([1,1:end-1],:) - x;
res(1,:) = -x(1,:);
res(end,:) = x(end-1,:);

function res = adjD3Dy(x)
res = x(:,[1,1:end-1],:) - x;
res(:,1,:) = -x(:,1,:);
res(:,end,:) = x(:,end-1,:);

function res = adjD3Dx(x)
res = x([1,1:end-1],:,:) - x;
res(1,:,:) = -x(1,:,:);
res(end,:,:) = x(end-1,:,:);

function res = adjD3Dz(x)
res = x(:,:,[1,1:end-1]) - x;
res(:,:,1) = -x(:,:,1);
res(:,:,end) = x(:,:,end-1);
