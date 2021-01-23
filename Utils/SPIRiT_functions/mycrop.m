function mat = mycrop(mat, cropsz, pos)
% Function res = MYCROP(mat, cropsize);
%
% Crop the input matrix to given size 
% From the first index (upper left corner in 2D)
%
% Stanislas Rapacchi, July 2011

% pos can be 'c'(centered) or 'ul'(upperleft) or 'dr' (downright)

if(~exist('pos','var'))
    pos='ul';
end

matsz=size(mat);
ndimIn=length(matsz);
ndimOut=length(cropsz);

if(ndimOut<ndimIn)
    Nmat=prod(matsz(ndimOut+1:ndimIn));
    cropsz=[cropsz(:)' matsz(ndimOut+1:ndimIn)];
%     omat=zeros([cropsz, Nmat]);
%     
%     for n=1:Nmat
%         imat=mat(
%         mat = mycrop(imat, cropsz, pos);
%     end
%     mat=omat;
else
    if(ndimOut>ndimIn)
        error('Too much crop dimensions');
        return;
    end
end

for i=1:length(cropsz)
    if(cropsz(i)<=1)
        cropsz(i)=matsz(i).*cropsz(i);
    end
end

bg=zeros(ndimOut,1);

switch pos
  case {'c'}
      for i=1:ndimOut
          bg(i)=1+floor((matsz(i)-cropsz(i))/2);
      end
  case {'dr'}
      for i=1:ndimOut
          bg(i)=1+(matsz(i)-cropsz(i));
      end
  case {'ul'}
      bg=ones(ndimOut,1);
  otherwise
      bg=ones(ndimOut,1);
end

%%
for dim=1:ndimOut
    mat=getindexindim(mat,dim,bg(dim):bg(dim)-1+cropsz(dim));
end
