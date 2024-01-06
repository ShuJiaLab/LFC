function A = getarea(imstackproj,rowc,colc)
[rows,cols] = ind2sub(size(imstackproj),find(imstackproj>0));
dists = sqrt((rows-rowc).^2 + (cols-colc).^2);

end