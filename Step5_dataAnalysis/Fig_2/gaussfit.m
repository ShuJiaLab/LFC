function [positions,sizes,mseErr] = gaussfit(Data2Fit,FittingMode,PixelSize)
%GAUSS Summary of this function goes here
%   Detailed explanation goes here
% vowel = vowel/max(vowel(:));
% Gauss2Mod = @(A,P1,P2,P3,XC,YC,ZC,k,x)...      % P is related to the size of the peak, greater->narrower
%     A*exp( -( 1./P1.*((x-XC).^2) + 1./P2.*((y-YC).^2) + 1./P3.*((z-ZC).^2)))+k;

GaussMod1D = @(BETA,x)...      % P is related to the size of the peak, greater->narrower
            BETA(1)*( exp( - 1./BETA(2).*((x(:,1)-BETA(3)).^2 )) )+BETA(4);
GaussMod2D = @(BETA,x)...      % P is related to the size of the peak, greater->narrower
            BETA(1)*(...
            exp( - 1./BETA(2).*((x(:,1)-BETA(3)).^2 )) .* ...
            exp( - 1./BETA(4).*((x(:,2)-BETA(5)).^2 )) ...
            )+BETA(6);
GaussMod3D = @(BETA,x)...      % P is related to the size of the peak, greater->narrower
            BETA(1)*(...
            exp( - 1./BETA(2).*((x(:,1)-BETA(3)).^2 )) .* ...
            exp( - 1./BETA(4).*((x(:,2)-BETA(5)).^2 )) .* ...
            exp( - 1./BETA(6).*((x(:,3)-BETA(7)).^2 ))   ...
            )+BETA(8);

Data2Fit = single(Data2Fit);
Data2Fit_proj = mean(Data2Fit,3);
Data2Fit_proj = Data2Fit_proj-min(Data2Fit_proj(:));
Data2Fit_proj = Data2Fit_proj/max(Data2Fit_proj(:));
[row,col,depth] = size(Data2Fit);
[maxv,maxloc] = max(max(Data2Fit,[],[1,2]));

a0 = maxv;
% figure(9999),imshow(Data2Fit_proj)
[rows,cols] = find(Data2Fit_proj>0.5*max(Data2Fit_proj(:)));
x0 = round(mean(cols,'all'));
y0 = round(mean(rows,'all'));
[~,maxi] = max(Data2Fit_proj(:));
[rowc,colc] = ind2sub(size(Data2Fit_proj),maxi);
if abs(x0-colc)>20 || abs(y0-rowc)>20
    x0 = colc;
    y0 = rowc;
end
rx0 = sqrt(length(cols));
ry0 = sqrt(length(rows));
rz0 = (rx0+ry0)/2;
% x0 = col/2;
% y0 = row/2;
z0 = maxloc;

if string(FittingMode) == "3D"
    StartP = [a0,(rx0^2)/log(16),x0,(ry0^2)/log(16),y0,(rz0^2)/log(16),z0,min(Data2Fit(:))];
%     disp(['Start Point: ',num2str(StartP(2:end))])
    [CoX,CoY,CoZ] = meshgrid(1:1:col,1:1:row,1:1:depth);
    X = [CoX(:),CoY(:),CoZ(:)];
    Y = single(Data2Fit(:));
    [beta,~,~,~,mseErr,~] = nlinfit(X,Y,GaussMod3D,StartP');
%     disp(['Fitted Params: ',num2str(beta(2:end)')])
    positions = [beta(3),beta(5),beta(7)];
    sizes = [sqrt(beta(2)*log(16))*PixelSize(1),sqrt(beta(4)*log(16))*PixelSize(2),sqrt(beta(6)*log(16))*PixelSize(3)];
elseif string(FittingMode) == "2D"
    StartP = [a0,(rx0^2)/log(16),x0,(ry0^2)/log(16),y0,min(Data2Fit(:))];
%     disp({'Start Point: ',StartP;'r: ',[rx0,ry0]})
    [CoX,CoY] = meshgrid(1:1:col,1:1:row);
    X = [CoX(:),CoY(:)];
    Y = single(Data2Fit(:));
    [beta,~,~,~,mseErr,~] = nlinfit(X,Y,GaussMod2D,StartP');
%     disp(['Fitted Params: ',num2str(beta')])
    positions = [beta(3),beta(5)];
    sizes = [sqrt(beta(2)*log(16))*PixelSize(1),sqrt(beta(4)*log(16))*PixelSize(2)];
elseif string(FittingMode) == "1D"
    vecind = find(Data2Fit_proj>0.75*max(Data2Fit_proj(:)));
    x0 = mean(vecind,'all');
    rx0 = length(vecind);
    StartP = [a0,(rx0^2)/log(16),x0,min(Data2Fit(:))];
%     disp(['Start Point: ',num2str(StartP(2:end))])
    CoX = 1:1:length(Data2Fit);
    X = [CoX(:)];
    Y = single(Data2Fit(:));
    [beta,~,~,~,mseErr,~] = nlinfit(X,Y,GaussMod1D,StartP');
%     disp(['Fitted Params: ',num2str(beta(2:end)')])
    positions = beta(3);
    sizes = sqrt(beta(2)*log(16))*PixelSize;
end
end

