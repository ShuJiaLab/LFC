function [imstack1out,imstack2out] = alignnumi(imstack1,imstack2)

%ALIGNNUMI Summary of this function goes here
%   Detailed explanation goes here

% imstack1 = magic(100);
% imstack2 = randn(100);

imstack1out = imstack1;
imstack2out = imstack2;

imstack1proj = mean(single(imstack1),3);
imstack2proj = mean(single(imstack2),3);

fig = uifigure("Name","Align nucleus to mitochondria","Position",[960, 540, 320, 330]);
ax = uiaxes('Parent',fig,'Units','pixels','Position', [0, 30, 300, 300]);
axis(ax,'tight');
imshowpair(imstack1proj,imstack2proj,"Parent",ax);
goupdown = uieditfield(fig,"numeric","Position",[30,25,100,24],"Value",0,"ValueChangedFcn",@moveimage);
goleftright = uieditfield(fig,"numeric","Position",[180,25,100,24],"Value",0,"ValueChangedFcn",@moveimage);

waitfor(fig);

    function moveimage(src,event)
        dy = goupdown.Value;
        dx = goleftright.Value;
        imstack1out = circshift(imstack1,[dy,dx]);
        imstack1proj = mean(single(imstack1out),3);
        imshowpair(imstack1proj,imstack2proj,"Parent",ax);drawnow;
    end
end