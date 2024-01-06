clear,clc
load("calc_Nusome_Dist/" + "NusomeDist_0.mat","NusomeDist");
NusomeDist_0 = NusomeDist(:,3);
load("calc_Nusome_Dist/" + "NusomeDist_30.mat","NusomeDist");
NusomeDist_30 = NusomeDist(:,3);
load("calc_Nusome_Dist/" + "NusomeDist_60.mat","NusomeDist");
NusomeDist_60 = NusomeDist(:,3);
load("calc_Nusome_Dist/" + "NusomeDist_120.mat","NusomeDist");
NusomeDist_120 = NusomeDist(:,3);
load("calc_Nusome_Dist/" + "NusomeDist_300.mat","NusomeDist");
NusomeDist_300 = NusomeDist(:,3);
clearvars NusomeDist
%%
dist_Nusome_0 = [];
for ii = 1:1:length(NusomeDist_0)
    dist_Nusome_0 = cat(1,dist_Nusome_0,NusomeDist_0{ii}');
end
dist_Nusome_30 = [];
for ii = 1:1:length(NusomeDist_30)
    dist_Nusome_30 = cat(1,dist_Nusome_30,NusomeDist_30{ii}');
end
dist_Nusome_60 = [];
for ii = 1:1:length(NusomeDist_60)
    dist_Nusome_60 = cat(1,dist_Nusome_60,NusomeDist_60{ii}');
end
dist_Nusome_120 = [];
for ii = 1:1:length(NusomeDist_120)
    dist_Nusome_120 = cat(1,dist_Nusome_120,NusomeDist_120{ii}');
end
dist_Nusome_300 = [];
for ii = 1:1:length(NusomeDist_300)
    dist_Nusome_300 = cat(1,dist_Nusome_300,NusomeDist_300{ii}');
end
%%
dist_cell_0 = [];
for ii = 1:1:length(NusomeDist_0)
    dist_cell_0 = cat(1,dist_cell_0,mean(NusomeDist_0{ii},"all"));
end
dist_cell_30 = [];
for ii = 1:1:length(NusomeDist_30)
    dist_cell_30 = cat(1,dist_cell_30,mean(NusomeDist_30{ii},"all"));
end
dist_cell_60 = [];
for ii = 1:1:length(NusomeDist_60)
    dist_cell_60 = cat(1,dist_cell_60,mean(NusomeDist_60{ii},"all"));
end
dist_cell_120 = [];
for ii = 1:1:length(NusomeDist_120)
    dist_cell_120 = cat(1,dist_cell_120,mean(NusomeDist_120{ii},"all"));
end
dist_cell_300 = [];
for ii = 1:1:length(NusomeDist_300)
    dist_cell_300 = cat(1,dist_cell_300,mean(NusomeDist_300{ii},"all"));
end