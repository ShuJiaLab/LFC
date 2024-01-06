filename = 'red_4000nm_2_6um_resized_clb.tif';
imstack = single(tiffreadVolume(filename));

imstack_proj1 = sum(imstack,3);
imstack_proj2 = squeeze(sum(imstack,1));

figure(1),imagesc(imstack_proj1);
[posx,posy] = ginput(1);
posx = round(posx);
posy = round(posy);

figure(1),imagesc(imstack_proj2);
[posz,pos0] = ginput(1);
posz = round(posz);
pos0 = round(pos0);
close("all")
%%
clc,
xprofile = imstack(posy,:,posz);xprofile = xprofile/max(xprofile);
yprofile = imstack(:,posx,posz);yprofile = yprofile/max(yprofile);

imstack_unified = imstack./max(imstack,[],[1,2]);
zprofile = squeeze(imstack_unified(posy,posx,:));
zprofile = zprofile/max(zprofile);
%%
figure(2),bx = bar(xprofile,'r');bx.FaceAlpha = 0.5;
hold on;
by = bar(yprofile,'g');by.FaceAlpha = 0.5;
figure(3),bz = bar(zprofile,'b');
%% +++++++++++++   Get Statistics   ++++++++++++++++++
% ====================================================
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
load("volstats_red.mat");volstats_red = volstats;
load("volstats_green.mat");volstats_green = volstats;
load("volstats_blue.mat");volstats_blue = volstats;
clearvars volstats
%%
figure(4),
scatter(cell2mat(volstats_red(:,4)),cell2mat(volstats_red(:,5)),'r','filled')
