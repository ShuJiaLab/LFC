clear,clc
%%
rawtif_folder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\XW20220618\z_multi_nucleus_acsn_ccut\";
nufolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\XW20220618\"+"PSFFLFint_20220630_Blue_Gly_10um-z_multi_nucleus_acsn_ccut";
prefix = "JK_minu60_";
filenames = dir(nufolder+"\*.tif");
filelist = {filenames(:).name}';
fileinds = split(filelist(:),[prefix,"(",")","_","iter50.tif"]);
filelist = [fileinds(:,1),num2cell(str2double(fileinds(:,[2,3,5]))),filelist];
filelist = sortrows(filelist);
filelist(:,1)={char(prefix)};
%% Open ImageJ
clc
warning off
addpath ('D:\FijiImageJ\scripts\');
addpath("Z:\Xuanwen\FLFMuf\ExpData\utilities\")
ImageJ;
%%
nufoldercropped = nufolder + "_cropped\";
nufolderlabeled = nufolder + "_labeled\";
mkdir(nufoldercropped)
mkdir(nufolderlabeled)
cropsizehalf = 200;
autothres = ["Default","Huang","Intermodes","IsoData","IJ_IsoData",...
             "Li","MaxEntropy","Mean","MinError","Minimum",...
             "Moments","Otsu","Percentile","RenyiEntropy","Shanbhag",...
             "Triangle","Yen"];
manualthres = 0;
%% Save cropped images
for ii =31:1:length(filelist)
    imstackname = nufolder +"\" + filelist{ii,5};
    imstack = single(tiffreadVolume(imstackname,"PixelRegion",{[1 1 inf],[1 1 inf],[1 1 inf]}));
    imstackproj = mean(imstack,3);
    [rows,cols] = find(imstackproj>(max(imstackproj(:))/2));
    rowc = round(mean(rows));
    colc = round(mean(cols));
    rowc = min([max([cropsizehalf,rowc]),size(imstack,1)-cropsizehalf]);
    colc = min([max([cropsizehalf,colc]),size(imstack,2)-cropsizehalf]);
    imstackcropped = imstack(rowc-cropsizehalf+1:rowc+cropsizehalf,colc-cropsizehalf+1:colc+cropsizehalf,:);
    imstackcroppedproj = mean(imstackcropped,3);
    imstackcropped = imstackcropped/max(imstackcropped(:));
    volwrite(imstackcropped*60000,nufoldercropped +"\" + filelist{ii,5})
    disp("Written " + filelist{ii,5})
end
%% 3D Nuclei Segmentation
clc
NusomeVol = [];
for ii = 19%31:1:length(filelist)
    fig = uifigure("WindowStyle","modal","Name","3D Segmentation","Position",[960, 540, 320, 450]);
    ax = uiaxes('Parent',fig,'Units','pixels','Position', [0, 130, 300, 300]);
    axis(ax,'tight');
    imshow(rescale(mean(single(imread(nufoldercropped +"\" + filelist{ii,5})),3),0,1),"Parent",ax);
    beginind = uieditfield(fig,"numeric","Position",[25,100,100,24],"Value",1);
    endind = uieditfield(fig,"numeric","Position",[175,100,100,24],"Value",101);
    atth = uicontrol(fig,"Style","popupmenu","Position",[25,60,100,24],"String",autothres);
    mnth = uieditfield(fig,"numeric","Position",[175,60,100,24],"Value",0);
    btn = uibutton(fig,'push','Position',[10, 20, 300, 24],...
               "Text","Segment",'ButtonPushedFcn',@(btn,event)...
               segButtonPushed(btn,nufoldercropped,nufolderlabeled,filelist{ii,5},atth,mnth,beginind,endind));
    waitfor(fig);
end
%% Save results
filenames = dir(nufolderlabeled+"\*.tif");
filelist = {filenames(:).name}';
fileinds = split(filelist(:),[prefix,"(",")","_","iter50.tif"]);
filelist = [fileinds(:,1),num2cell(str2double(fileinds(:,[2,3,5]))),filelist];
filelist = sortrows(filelist);
filelist(:,1)={char(prefix)};
%%
NusomeVol = [];
for ii =1:1:length(filelist) 
    if ~exist(nufolderlabeled+"\"+filelist{ii,5},"file")
        return
    else
        imstacklabeled = tiffreadVolume(nufolderlabeled+"\"+filelist{ii,5});
        numofnu = max(imstacklabeled(:));
        vols = zeros(1,numofnu);
        for jj = 1:1:numofnu
            vols(jj) = sum(double(imstacklabeled==jj),"all");
        end
        vols = double(vols)*0.065*0.065*0.1;
        vols(vols<pi/6*0.9^3)=[];
        NusomeVol = cat(1,NusomeVol,{filelist{ii,5},numofnu,vols,(vols*6/pi).^(1/3)});
        disp(ii+" ["+filelist{ii,5}+"] "+numofnu+": ["+num2str((vols*6/pi).^(1/3))+"] * Saved! *")
    end

%     rawim = filelist{ii,5};
%     rawim = single(imread(rawtif_folder + rawim(1:end-10) + ".tif"));
%     rawim = rawim(476:476+500,260:260+500);
%     subplot(121),imshow(rawim/max(rawim(:))),
%     title("["+filelist{ii,2}+" "+filelist{ii,3}+" "+filelist{ii,4}+"]")
%     subplot(122),scatter(vols*0.0 + double(numofnu),(vols*6/pi).^(1/3),'filled'),hold on
%     drawnow
end
% hold off

%% utilities
function segButtonPushed(~,nufoldercropped,nufolderlabeled,filename,atth,mnth,beginind,endind)
    autothres = ["Default","Huang","Intermodes","IsoData","IJ_IsoData",...
             "Li","MaxEntropy","Mean","MinError","Minimum",...
             "Moments","Otsu","Percentile","RenyiEntropy","Shanbhag",...
             "Triangle","Yen"];
    disp(autothres(atth.Value) + " " + mnth.Value)
    manualthres = mnth.Value; 
    ij.IJ.open(nufoldercropped+"\"+filename);
    ij.IJ.selectWindow(filename)
    ij.IJ.run("Slice Keeper", "first="+beginind.Value+" last="+endind.Value+" increment=1");
    ij.IJ.selectWindow(filename)
    ij.IJ.run("Close","");
    ij.IJ.selectWindow(filename+" kept stack")
    ij.IJ.run("3D Nuclei Segmentation",...
        "auto_threshold="+autothres(atth.Value)+" manual="+manualthres+" separate_nuclei");
    ij.IJ.selectWindow("merge");
    ij.IJ.saveAs("Tiff",nufolderlabeled+"\"+filename);
    ij.IJ.selectWindow(filename+" kept stack")
    ij.IJ.run("Close","");

    imstacklabeled = tiffreadVolume(nufolderlabeled+"\"+filename);
    numofnu = max(imstacklabeled(:));
    vols = zeros(1,numofnu);
    for jj = 1:1:numofnu
        vols(jj) = sum(double(imstacklabeled==jj),"all");
    end
    disp("["+filename+"] "+numofnu+": Vol(px) ["+vols+"]")
    vols = double(vols)*0.065*0.065*0.1;
    vols(vols<pi/6*0.9^3)=[];
    disp("["+filename+"] "+numofnu+": Diamater ["+num2str((vols*6/pi).^(1/3))+"]")
end
