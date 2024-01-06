function [tomove,z0] = zprof_select(imstack1,imstack2)
%ZPROF_SELECT Summary of this function goes here
%   Detailed explanation goes here
zprofile_max_1 = rescale(squeeze(max(imstack1,[],[1,2])),0,1);
zprofile_sum_1 = rescale(squeeze(sum(imstack1,[1,2])),0,1);

zprofile_max_2 = rescale(squeeze(max(imstack2,[],[1,2])),0,1);
zprofile_sum_2 = rescale(squeeze(sum(imstack2,[1,2])),0,1);

tomove = 0;
z0 = 30;
fig = uifigure("Name","Align nucleus to mitochondria","Position",[960, 540, 640, 380]);
ax1 = uiaxes('Parent',fig,'Units','pixels','Position', [0, 50, 300, 300]);axis(ax1,'tight');
ax2 = uiaxes('Parent',fig,'Units','pixels','Position', [320, 50, 300, 300]);axis(ax2,'tight');
plot(1:1:length(zprofile_max_1),zprofile_max_1,1:1:length(zprofile_sum_1),zprofile_sum_1,"Parent",ax1);
plot(1:1:length(zprofile_max_2),zprofile_max_2,1:1:length(zprofile_sum_2),zprofile_sum_2,"Parent",ax2);
focalplane = uieditfield(fig,"numeric","Position",[30,10,100,24],"Value",30,"ValueChangedFcn",@selectz);
tomovebox = uicheckbox(fig,"Text","Move to trash?","Position",[180,10,100,24],"Value",0,"ValueChangedFcn",@movechanged);

waitfor(fig);

    function movechanged(src,event)
        tomove = tomovebox.Value;
    end

    function selectz(src,event)
        z0 = focalplane.Value;
    end

end

