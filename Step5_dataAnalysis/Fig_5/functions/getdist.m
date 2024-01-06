function dist = getdist(imstack,imstackvolrep,positive,rowc,colc,pagec)
[rows,cols,pages] = ind2sub(size(imstack),find(imstack>positive));
dists = sqrt((rows-rowc).^2 + (cols-colc).^2 + (pages-pagec).^2);
dists = (dists/imstackvolrep);
% dists = dists/max(dists(:));
% disp(length(dists))
[counts,edges] = histcounts(dists,50);
% Q = uniformspherepdf(edges,imstackvolrep);
% Q = counts.*0 + 1/length(counts);
% KL = KLdiv(counts,Q);
edges = (edges(1:end-1) + edges(2:end))/2;
% figure(3),plot(edges,counts)
f = fit(edges',counts','gauss1');
dist = f.b1;
end