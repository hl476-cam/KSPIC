function col = dscatter_col(X,Y)
lambda = [];
nbins = [];
plottype = 'scatter';
contourFlag = false;
% msize = 10;
% marker = 's';


minx = min(X,[],1);
maxx = max(X,[],1);
miny = min(Y,[],1);
maxy = max(Y,[],1);

if isempty(nbins)
    nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
end

if isempty(lambda)
    lambda = 20;
end

edges1 = linspace(minx, maxx, nbins(1)+1);
ctrs1 = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(miny, maxy, nbins(2)+1);
ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];

[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
[dum,bin(:,2)] = histc(X,edges1);
[dum,bin(:,1)] = histc(Y,edges2);
H = accumarray(bin,1,nbins([2 1])) ./ n;
G = smooth1D(H,nbins(2)/lambda);
F = smooth1D(G',nbins(1)/lambda)';
% = filter2D(H,lambda);

okTypes = {'surf','mesh','contour','image','scatter'};
k = strmatch(lower(plottype), okTypes); %#ok
if isempty(k)
    error('dscatter:UnknownPlotType',...
        'Unknown plot type: %s.',plottype);
elseif length(k)>1
    error('dscatter:AmbiguousPlotType',...
        'Ambiguous plot type: %s.',plottype);
else
    switch(k)
        
        case 1 %'surf'
            h = surf(ctrs1,ctrs2,F,'edgealpha',0);
        case 2 % 'mesh'
            h = mesh(ctrs1,ctrs2,F);
        case 3 %'contour'
            [dummy, h] =contour(ctrs1,ctrs2,F);
        case 4 %'image'
            nc = 256;
            F = F./max(F(:));
            colormap(repmat(linspace(1,0,nc)',1,3));
            h =image(ctrs1,ctrs2,floor(nc.*F) + 1);
        case 5 %'scatter'
            F = F./max(F(:));
            ind = sub2ind(size(F),bin(:,1),bin(:,2));
            col = F(ind);

%             figure;
%             h = scatter(X,Y,msize,col.^SCALE,marker,'filled');colormap(jet);caxis([0 BAR]);axis square;
            
    end
    
end



%%%% This method is quicker for symmetric data.
% function Z = filter2D(Y,bw)
% z = -1:(1/bw):1;
% k = .75 * (1 - z.^2);
% k = k ./ sum(k);
% Z = filter2(k'*k,Y);

function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
end

end

