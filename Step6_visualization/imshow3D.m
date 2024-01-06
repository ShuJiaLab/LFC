function camprops = imshow3D(imstack,im3dsavingpath,camadj,varargin)
% Show 3D stacks using boxed volume rendering. The function convert the 3D 
% volume to a volume rendered view with a certain camera orientation using
% Maximum Intensity Projection algorithm.
%
%# 3D shown by setting and return the camera orientation
% camprops = imshow3D(imstack,im3dsavingpath,true)
%
%# 3D shown using existing camera orientation
% camprops = imshow3D(imstack,im3dsavingpath,false,camprops0)
% -----------------------------------------------
% Input parameters
%  - imstack: input 3D volume
%  - im3dsavingpath: target location to save the rendered view in a single
%  2D RGB image
%  - camadj: bool value to 
%  - varargin{1}: set if to show 3D rendered view using existing camera
%  orientation. The camera orientation is kept in the format 
%  [{"CameraPosition",campos},{"CameraTarget",camtarget},{"CameraUpVector",
%  camvec,"CameraZoom",camzoom}].
% -----------------------------------------------
% Output parameters
%  - camprops: output camera orientation.
% -----------------------------------------------
% Copyright 2022 Xuanwen Hua
% x.hua@gatech.edu

import java.awt.Robot;
mouse = Robot;
camprops = [];

fig2 = uifigure("Name","Volume-Only View");
fig2.Position = [500,400,500,500];
viewer2 = viewer3d(fig2);
viewer2.BackgroundGradient = "off";
viewer2.BackgroundColor = "black";
viewer2.Box = "on";
viewer2.ScaleBar = "on";
viewer2.ScaleBarUnits = "voxels = "+500*0.065+" μm";
viewer2.Lighting = "off";
vol2 = volshow(imstack,Parent=viewer2);
vol2.RenderingStyle = "MaximumIntensityProjection";

if camadj
    fig = uifigure("Name","Volume Orientation Pre-setting"+" ("+num2str(size(imstack))+")");
    fig.Position = [1000,400,500,500];
    grid = uigridlayout(fig);
    grid.RowHeight = {'1x',20,20,20,20};
    grid.ColumnWidth = {'1x'};
    viewer = viewer3d(grid);
    viewer.BackgroundGradient = "off";
    viewer.BackgroundColor = "black";
    viewer.Box = "on";
    viewer.ScaleBar = "on";
    viewer.ScaleBarUnits = "voxels = "+500*0.065+" μm";
    viewer.Lighting = "off";
    vol = volshow(imstack,Parent=viewer);
    vol.RenderingStyle = "MaximumIntensityProjection";
    % vol.GradientOpacityValue = 0.1;
    % --------------------------------------------------
    campos = uilabel(grid,"Text","Camera Position: -- -- --");
    camtarget = uilabel(grid,"Text","CameraTarget: -- -- --");
    camupvec = uilabel(grid,"Text","CameraUpVector: -- -- --");
    camzoom = uilabel(grid,"Text","CameraZoom: -- -- --");
    camlabels = [campos,camtarget,camupvec,camzoom];
    camtypes = ["CameraPosition","CameraTarget","CameraUpVector","CameraZoom"];
    addlistener(viewer,"CameraMoving",@(scr,evt)printcampos(viewer,camlabels,camtypes,viewer2));
    addlistener(viewer,"CameraMoved",@(scr,evt)printcampos(viewer,camlabels,camtypes,viewer2));
    waitfor(fig);
    for ii = 1:1:length(camtypes)
        camtype = camtypes(ii);
        camprops = cat(1,camprops,{camtype,get(viewer2,camtype)});
    end
    mouse.mouseMove(900, 300);
    exportapp(fig2,im3dsavingpath);
    pause(0.1);
    delete(fig2);
    mouse.mouseMove(1100, 100);
else
    camprops = varargin{1};
    for ii = 1:1:length(camprops)
        set(viewer2,camprops{ii,1},camprops{ii,2});
    end    
    pause(2);
    mouse.mouseMove(900, 300);
    pause(1);
    mouse.mouseMove(905, 305);
    exportapp(fig2,im3dsavingpath);
    pause(1);
    delete(fig2);
end
end
%%
function printcampos(scr,labelctrls,camtypes,viewer2)
for ii = 1:1:length(labelctrls)
    labelctrl = labelctrls(ii);
    camtype = camtypes(ii);
    labelctrl.Text = camtype + ": "+ num2str(get(scr,camtype));
    labelctrl.FontSize = 12;
    set(viewer2,camtype,get(scr,camtype));
end
end