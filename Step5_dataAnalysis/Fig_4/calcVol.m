clear,clc
addpath('Z:\Xuanwen\FLFMuf\ExpData\utilities\')

mbfolder = './PSFFLFint_Sim65nm_20220206_Red_refine_gly_10um_1024-XW20220215_rawtif_selected_mb_acsn_bkgsub_ccut/';
nufolder = './PSFFLFint_Sim65nm_20220206_Blue_gly_10um_1024-XW20220215_rawtif_selected_nu_acsn_bkgsub_ccut/';
mbdir = dir([mbfolder,'*.tif']);
nudir = dir([nufolder,'*.tif']);
%%
mbnu_vrate = [];
H = waitbar(0,'Calculating...');
for ii = 1:1:length(mbdir)
mbstack = tiffreadVolume([mbdir(ii).folder,'\',mbdir(ii).name]);
mbstack = double(mbstack(611:1710,591:1690,:));
mbstack = mbstack/max(mbstack(:));
% figure(1),imshow(mbstack(:,:,51))
mbV = calcVol_mb(mbstack,21,81);

nustack = tiffreadVolume([nudir(ii).folder,'\',nudir(ii).name]);
nustack = double(nustack(611:1710,591:1690,:));
nustack = nustack/max(nustack(:));
% figure(2),imshow(nustack(:,:,51))
nuV = calcVol_nu(nustack,21,81);

mbnu_vrate = cat(1,mbnu_vrate,[mbV,nuV,nuV/mbV]);
waitbar(ii/length(mbdir),H,['# ',num2str(ii),'/',num2str(length(mbdir)),', [',num2str([mbV,nuV,nuV/mbV]),']'])
end
disp('VOL calculation done!')
delete(H)
%%
% [X,Y,Z] = meshgrid(1:1:size(imstackcut,1),1:1:size(imstackcut,2),1:1:size(imstackcut,3));
% xc = sum(X.*imstackcut,'all')/sum(imstackcut,'all');
% yc = sum(Y.*imstackcut,'all')/sum(imstackcut,'all');
% zc = sum(Z.*imstackcut,'all')/sum(imstackcut,'all');
%%
imstackcut = mbstack;
% imstackcut = imstackcut/max(imstackcut(:));
clc,
% figure(6),
imV = 0;
H = waitbar(0,'Calculating volume...');
for depth = 21%21:1:81%size(imstackcut,3)
    imlayer = imstackcut(:,:,depth);
    [~,threshold] = edge(imlayer,'sobel');
    fudgeFactor = 5;
    BWs = edge(imlayer,'sobel',threshold * fudgeFactor);
    subplot(231),imshow(BWs),title('Binary Gradient Mask')
    strelsize = 5;
    se135 = strel('line',strelsize,135);
    se90 = strel('line',strelsize,90);
    se45 = strel('line',strelsize,45);
    se0 = strel('line',strelsize,0);
    BWsdil = imdilate(BWs,[se135 se90 se45 se0]);
    subplot(232),imshow(BWsdil),title('Dilated Gradient Mask')
%     BWdfill = imfill(BWsdil,'holes');
    BWdfill = imclose(BWsdil,strel('disk',30));
    subplot(233),imshow(BWdfill),title('Binary Image with Filled Holes')
    BWnobord = imclearborder(BWdfill,4);
    subplot(234),imshow(BWnobord),title('Cleared Border Image')
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    subplot(235),imshow(BWfinal),title('Segmented Image');
    BWoutline = bwperim(BWfinal);
    Segout = imstackcut(:,:,depth)/max(max(imstackcut(:,:,depth)));
    Segout(BWoutline) = 1;
    subplot(236),imshow(Segout*255,hot),title('Outlined Original Image')
    imV = imV + sum(BWfinal,'all');
    waitbar((depth-21)/61,H,'Calculating...')
end
delete(H)
%%
mbnu_vrate_real = mbnu_vrate;
mbnu_vrate_real(:,1:2) = (mbnu_vrate_real(:,1:2)*(0.1*0.065^2));