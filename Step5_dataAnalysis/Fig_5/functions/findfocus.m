function [rowc, colc,pagec] = findfocus(imstack)
energygradientx = (imstack(2:end,:,:)-imstack(1:end-1,:,:)).^2;
energygradienty = (imstack(:,2:end,:)-imstack(:,1:end-1,:)).^2;
pagec = sum(energygradientx,[1,2]) + sum(energygradienty,[1,2]);
pagec  = squeeze(pagec);
% figure(101),plot(pagec);
rowc = 1;
colc = 1;
end
