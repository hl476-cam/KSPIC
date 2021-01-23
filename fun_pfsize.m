function pfSize = fun_pfsize(pattern)

% A simple but not very elegant way to detect the partial Fourier size

pattern=double(pattern);

pfSize=size(pattern);

for iRow=round(size(pattern,1)/2):size(pattern,1)-3
    row=[pattern(iRow,:);pattern(iRow+1,:);pattern(iRow+2,:);pattern(iRow+3,:)]; % when the next 4 lines are all empty
    if norm(row,2)==0
        pfSize(1)=iRow-1;
        break;
    end  
end

for iColumn=round(size(pattern,2)/2):size(pattern,2)-3
    column=[pattern(:,iColumn);pattern(:,iColumn+1);pattern(:,iColumn+2);pattern(:,iColumn+3)]; % when the next 4 lines are all empty
    if norm(column,2)==0
        pfSize(2)=iColumn-1;
        break;
    end  
end

end

