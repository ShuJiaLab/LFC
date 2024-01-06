function cpstack = cropstack(imstack)
cropsizehalf = 200;
imstack = single(imstack);
imstackproj = mean(imstack,3);
[rows,cols] = find(imstackproj>(max(imstackproj(:))*0.5));
rowc = round(mean(rows));
colc = round(mean(cols));
rowc = min([max([cropsizehalf,rowc]),size(imstack,1)-cropsizehalf]);
colc = min([max([cropsizehalf,colc]),size(imstack,2)-cropsizehalf]);
cpstack = imstack(rowc-cropsizehalf+1:rowc+cropsizehalf,colc-cropsizehalf+1:colc+cropsizehalf,:);
end