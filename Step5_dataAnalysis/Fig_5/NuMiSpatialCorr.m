clear,clc
addpath("Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\functions");
addpath("Z:\Xuanwen\FLFMuf\ExpData\utilities\")
params = ["XW20220609","XW20220616","XW20220618","XW20220619","XW20220620";
    "JK_minu","JK_minu30_","JK_minu60_","JK_minu120_","JK_minu300_"];

groupind = 5;
groupname = params(1,groupind) ;
prefix = params(2,groupind);
nufolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\"+groupname+"\PSFFLFint_20220630_Blue_Gly_10um-z_multi_nucleus_acsn_ccut";
mifolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\"+groupname+"\PSFFLFint_20220630_Red_Gly_10um-z_multi_mito_acsn_ccut";

nufilenames = dir(nufolder+"\*.tif");
nufilelist = {nufilenames(:).name}';
nufileinds = split(nufilelist(:),[prefix,"(",")","_","iter50.tif"]);
nufilelist = [nufileinds(:,1),num2cell(str2double(nufileinds(:,[2,3,5]))),nufilelist];
nufilelist = sortrows(nufilelist);

mifilenames = dir(mifolder+"\*.tif");
mifilelist = {mifilenames(:).name}';
mifileinds = split(mifilelist(:),[prefix,"(",")","_","iter80.tif"]);
mifilelist = [mifileinds(:,1),num2cell(str2double(mifileinds(:,[2,3,5]))),mifilelist];
mifilelist = sortrows(mifilelist);

numifilelist = [nufilelist(:,5),mifilelist(:,5)];
removelist = [];

% make directories
nusavefolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\"+groupname+"\000_mi_nu";
misavefolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\"+groupname+"\000_mi_all";
warning off;
mkdir(nusavefolder);mkdir(misavefolder)
%% New
clc;
dutycycle = [];
for ii = 19%1:100
    nustack = cropstack(tiffreadVolume(nufolder + "\"+numifilelist{ii,1},"PixelRegion",{[1 1 inf],[1 1 inf],[11 1 90]}));
    mistack = cropstack(tiffreadVolume(mifolder + "\"+numifilelist{ii,2},"PixelRegion",{[1 1 inf],[1 1 inf],[11 1 90]}));
    disp(ii+" | "+numifilelist{ii,1}+", "+numifilelist{ii,2})
    [nustackaligned,mistackaligned,rowc,colc,pagec,tomove,...
        nustacksrf,mistacksrf,mistack_nu,mistack_all] = NuMiCore(nustack,mistack);
    if tomove==1
        removelist = cat(1,removelist,[numifilelist(ii,1),numifilelist(ii,2)]);
        disp("[~] " + numifilelist{ii,1}+", "+numifilelist{ii,2} + " [excluded]")
    else
%         volwrite(mistack_nu,nusavefolder+"\"+numifilelist{ii,2})
%         volwrite(mistack_all,misavefolder+"\"+numifilelist{ii,2})
        disp("Current items: (Nu) "+length(dir(nusavefolder+"\*.tif"))+" | (Mi) "+length(dir(misavefolder+"\*.tif")))
        dutycycle = cat(1,dutycycle,sum(double(nustacksrf).*double(mistacksrf),"all")/sum(double(nustacksrf),"all"));
        disp("Duty Cycle: "+dutycycle);
    end
end
%% ===================================================================================================
% ====================================================================================================
clc;voldc = [];
nusavefilenames = dir(nusavefolder+"\*.tif");
nusavefilelist = {nusavefilenames(:).name}';
nusavefileinds = split(nusavefilelist(:),[prefix,"(",")","_","iter80.tif"]);
nusavefilelist = [nusavefileinds(:,1),num2cell(str2double(nusavefileinds(:,[2,3,5]))),nusavefilelist];
nusavefilelist = sortrows(nusavefilelist);

for ii = 1:1:length(nusavefilelist)
    mistack_nu = tiffreadVolume(nusavefolder+"\"+ nusavefilelist{ii,5});
    mistack_all = tiffreadVolume(misavefolder+"\"+ nusavefilelist{ii,5});
    voldc = cat(1,voldc,[sum(double(mistack_nu>0),"all"),sum(double(mistack_all>0),"all"),...
        sum(double(mistack_nu>0).*double(mistack_all>0),"all"),...
        sum(double((mistack_nu+mistack_all)>0),"all")]);
    disp("Duty Cycle: "+voldc(ii,:));
end