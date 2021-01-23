function mat=getindexindim(mat,dim,index)
% FUNCTION mat=GETINDEXINDIM(mat,dim,index)
% to apply a vector to number (like max, mean, etc) operation 
% 'operation' on a single dimension of matrix in_mat
% Output is the matrix with chosen dimension equal to 1
% EXAMPLES of usage:
% load wmri
% myoperation=@(x)max(x)/sum(x);
% outmat=concatenate_operation(X,3,myoperation);
%
% Perform a root-sum of square on a 4D matrix on the 4th dimension:
% im3D = sqrt(concatenate_operation(mat4D.^2,4,@sum));
%
% Stanislas Rapacchi
% July 2011

inmatsz = size(mat);

if(dim>length(inmatsz) || dim<1)
    error('incorrect chosen dimension');
end


% put the chosen dim as the first dim
if(dim > 1)
    permutesz = [dim (1:dim-1) (dim+1:size(inmatsz,2))];
    mat = permute(mat,permutesz);
else
    permutesz = 1:length(inmatsz);
    mat = mat;
end

% get the indexes
tempsz=inmatsz(permutesz);
tempsz(1)=length(index);

mat = reshape(mat(index,:),tempsz);



%permute back
if(dim > 1 )
    if(dim<length(inmatsz))
        invpermutesz = [(2:dim) 1 (dim+1:length(inmatsz))];
    else
        invpermutesz = [(2:length(inmatsz)) 1 ];
    end
    mat = permute(mat,invpermutesz);
else
    mat = mat;
end

return;
