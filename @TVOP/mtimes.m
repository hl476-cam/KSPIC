function res = mtimes(a,b)

minus1D = a.minus1D;

%wTV
bsz=size(b);
TVtype=a.TVtype;
aTV=a.weight;

if a.adjoint
	res = adjD(b,minus1D);

    %wTV
     %res = iwTV(b);

else
	res = D(b,minus1D);

    %wTV
     %res = zeros([bsz,TVtype/2]);
     %for x=1:bsz(3)
      %   res(:,:,x,:)=weightedTV2D(b(:,:,x),TVtype,aTV);
     %end

end




    
