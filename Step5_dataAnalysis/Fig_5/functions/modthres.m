function thres = modthres(thres0,zsize,FWHM)
%MODTHRES Summary of this function goes here
%   Detailed explanation goes here
X = 1:1:zsize;
miu = round(zsize/2);
sigma = ((FWHM/2)^2)/(log(2));
modTF = exp((-(X-miu).^2)/sigma);
thres = 1 - (1-thres0).* modTF;
thres = reshape(thres,1,1,[]);
end

