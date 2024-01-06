clear,clc
load("NusomeVol_0.mat","NusomeVol");
NusomeVol_0 = NusomeVol(:,3);
load("NusomeVol_30.mat","NusomeVol");
NusomeVol_30 = NusomeVol(:,3);
load("NusomeVol_60.mat","NusomeVol");
NusomeVol_60 = NusomeVol(:,3);
load("NusomeVol_120.mat","NusomeVol");
NusomeVol_120 = NusomeVol(:,3);
load("NusomeVol_300.mat","NusomeVol");
NusomeVol_300 = NusomeVol(:,3);

%%
broken = [];
unbroken = [];
NusomeVol = NusomeVol_300;
for ii = 1:1:length(NusomeVol)
    if length(NusomeVol{ii})==1
        unbroken = cat(2,unbroken,NusomeVol{ii});
    else
        broken = cat(2,broken,NusomeVol{ii});
    end
end
%%
stats = cat(1,stats,[300,length(unbroken),mean(unbroken),std(unbroken),length(broken),mean(broken),std(broken)]);
%%
clc
szfactor = 10;
figure(5);
for jj = 1:1:size(stats,1)
    scatter(stats(jj,1),stats(jj,3),stats(jj,2)*szfactor,"blue","filled");hold on;
    scatter(stats(jj,1),stats(jj,6),stats(jj,5)*szfactor,"red","filled");hold on;
end
hold off
%% get nucleosome numbers
NumNusome_300 = [];
for ii = 1:1:length(NusomeVol_300)
    NumNusome_300 = cat(1,NumNusome_300,length(NusomeVol_300{ii}));
end