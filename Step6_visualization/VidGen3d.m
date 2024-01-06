clear,clc
% -----------------------------------------------------------
% Get 3D stack lists
stackname = '.\JurkatSTS\JK_minu30_00001(5)_250_291.tif';
reconstackfolders = dir("C:\Users\Xuanwen\Desktop\videos\PSF*");
iternum = 50;
stacknamecut = split(stackname,["_","."]);
startind = str2double(stacknamecut{end-2});
endind = str2double(stacknamecut{end-1});
reconstacks_all = [];
for ii = 1:1:length(reconstackfolders)
    reconstacks_all = cat(1,reconstacks_all,dir([reconstackfolders(ii).folder,'\',reconstackfolders(ii).name,'\*.tif']));
end

stacklist = [];
for ii = startind:1:endind
    stacknamecut = split(stackname,["\",")"]);
    reconstaskname = [stacknamecut{end-1},')_',num2str(ii),'iter',num2str(iternum),'.tif'];
    reconstackpath = reconstacks_all(string({reconstacks_all(:).name}')==string(reconstaskname));
    reconstackfolder = reconstackpath.folder;
    reconstackfoldercut = split(reconstackfolder,'_');
    if ismember({'Red'},reconstackfoldercut)
        stacklist = cat(1,stacklist,[string(reconstackpath.folder),string(reconstackpath.name),"Red"]);
    elseif ismember({'Blue'},reconstackfoldercut)
        stacklist = cat(1,stacklist,[string(reconstackpath.folder),string(reconstackpath.name),"Blue"]);
    end
end
%% Show 3D
stacknamecut = split(stackname,["\","."]);
savefolder = fileparts(stacklist(1,1))+"\V_"+stacknamecut(end-1)+"\";
mkdir(savefolder);
camprops = [];
for ii = 1:1:length(stacklist)
    imstack = tiffreadVolume(stacklist(ii,1)+"\"+stacklist(ii,2));
    [rows,cols,deps] = size(imstack);
    imstack = double(imresize3(imstack,round([rows,cols,deps*100/65])));
    imstack = rescale(imstack,0,255);
    if stacklist(ii,3)=="Red"
        imstack = cat(4,imstack,imstack*0,imstack);
    elseif stacklist(ii,3)=="Green"
        imstack = cat(4,imstack*0,imstack,imstack*0);
    elseif stacklist(ii,3)=="Blue"
        imstack = cat(4,imstack*0,imstack,imstack);
    end
    imstack = uint8(imstack);
    camprops = imshow3D(imstack,savefolder+"V_"+stacklist(ii,2),ii==1,camprops);
end

%%
imstack = tiffreadVolume('testvol.tif');
rows = size(imstack,1);
cols = size(imstack,2);
deps = size(imstack,3);
imstack = cat(4,imresize3(imstack(:,:,:,1),round([rows,cols,deps*100/65])),...
                imresize3(imstack(:,:,:,2),round([rows,cols,deps*100/65])),...
                imresize3(imstack(:,:,:,3),round([rows,cols,deps*100/65])));
disp(size(imstack))
% imshow(squeeze(imstack(:,:,51,:)))
fig = uifigure("Name","Volume Orientation Pre-setting");
fig2 = uifigure("Name","Volume-Only View");
fig.CloseRequestFcn = @(scr,evt)myclosefcn(fig,fig2);
grid = uigridlayout(fig);
grid.RowHeight = {'1x',20,20,20,20};
grid.ColumnWidth = {'1x'};
viewer = viewer3d(grid);
viewer2 = viewer3d(fig2);
% --------------------------------------------------
viewer.BackgroundGradient = "off";
viewer.BackgroundColor = "black";
viewer.Box = "on";
viewer.ScaleBar = "on";
viewer.ScaleBarUnits = "voxels = 6.5 μm";
viewer.Lighting = "off";
viewer2.BackgroundGradient = "off";
viewer2.BackgroundColor = "black";
viewer2.Box = "on";
viewer2.ScaleBar = "on";
viewer2.ScaleBarUnits = "voxels = 6.5 μm";
viewer2.Lighting = "off";

vol = volshow(imstack,Parent=viewer);
vol.RenderingStyle = "MaximumIntensityProjection";
% vol.GradientOpacityValue = 0.1;
vol2 = volshow(imstack,Parent=viewer2);
vol2.RenderingStyle = "MaximumIntensityProjection";
% --------------------------------------------------
campos = uilabel(grid,"Text","Camera Position: -- -- --");
camtarget = uilabel(grid,"Text","CameraTarget: -- -- --");
camupvec = uilabel(grid,"Text","CameraUpVector: -- -- --");
camzoom = uilabel(grid,"Text","CameraZoom: -- -- --");
camlabels = [campos,camtarget,camupvec,camzoom];
camtypes = ["CameraPosition","CameraTarget","CameraUpVector","CameraZoom"];
addlistener(viewer,"CameraMoving",@(scr,evt)printcampos(viewer,camlabels,camtypes,viewer2));
addlistener(viewer,"CameraMoved",@(scr,evt)printcampos(viewer,camlabels,camtypes,viewer2));

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

function myclosefcn(scr,fig2)
try
exportapp(fig2,"test.png")
delete(fig2)
catch
end
delete(scr)
end