function [row,col,page] = getcent(imstack,thres)
% [~,threshold] = edge(imstack(:,:,80),'sobel');
[rows,cols,pages] = ind2sub(size(imstack),find(imstack>(thres*max(imstack,[],[1,2]))));
row = mean(rows,"all");
col = mean(cols,"all");
page = mean(pages,"all");
end