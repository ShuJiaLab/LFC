clear,clc
load("res.mat","res"); 
radii_relative = [1-res(:,2)./res(:,1),1-res(:,3)./res(:,1),res(:,4)];
% foldername = "./PSFFLFint_20220515_Blue_Gly_10um-rawtif_selected_focused_offsetsub_acsn_ccut/";
% filelist = dir(foldername+"/*.tif");
%%
cluster_num = 3;
clc,[idx,C] = kmeans(radii_relative(:,1:3),cluster_num,"MaxIter",100);
%
for ii = 1:100
    [idx,C] = kmeans(res(:,1:3),cluster_num,"MaxIter",1000,"Start",C);
end
figure(4444),gscatter(1-res(:,2)./res(:,1),1-res(:,3)./res(:,1),idx,'bgmc')
hold on
plot(1-C(:,2)./C(:,1),1-C(:,3)./C(:,1),'kx'),drawnow;
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster Centroid')

% e1se2s = [];
% e1se2l = [];
% e1le2s = [];
% e1le2l = [];
% for ii = 1:1:100
%     mbstack = single(tiffreadVolume(fullfile(filelist(ii).folder,filelist(ii).name)));
%     mbstackproj = rescale(max(mbstack,[],3),0,1);
%     switch idx(ii)
%         case 1
%             e1se2s = cat(2,e1se2s,mbstackproj);
%         case 2
%             e1se2l = cat(2,e1se2l,mbstackproj);
%         case 3
%             e1le2s = cat(2,e1le2s,mbstackproj);
% %         case 4
% %             e1le2l = cat(2,e1le2l,mbstackproj);
%     end
%     disp("#"+ii+" "+filelist(ii).name)
% end
% figure(4445),
% subplot(411),imshow(e1se2s*255,hot(256)),title("[Group #1]");
% subplot(412),imshow(e1se2l*255,hot(256)),title("[Group #2]");
% subplot(413),imshow(e1le2s*255,hot(256)),title("[Group #3]");
% subplot(414),imshow(e1le2l*255,hot(256)),title("[Group #4]");

%%
foldername = "./PSFFLFint_20220515_Blue_Gly_10um-rawtif_selected_focused_offsetsub_acsn_ccut/";
filelist = dir(foldername+"/*.tif");
load("res_k3.mat","res"); 
% res = [a,b,c,v,idx]; radii_relative = [e1,e2,V];
radii_relative = [1-res(:,2)./res(:,1),1-res(:,3)./res(:,1),res(:,4)];
colidx = repmat([1,0,1],length(radii_relative),1);
colidx(res(:,5)==2,:) = colidx(res(:,5)==2,:)*0 + [0.1,0.75,0.1];
colidx(res(:,5)==3,:) = colidx(res(:,5)==3,:)*0 + [0,0.5,1];
figure(4444),scatter3(res(:,1),res(:,2),res(:,3),30,colidx,"filled")

res_g1 = res(res(:,5)==1,:);
res_g2 = res(res(:,5)==2,:);
res_g3 = res(res(:,5)==3,:);
% radii_relative_gi = [e1,e2,V,ind];
radii_relative_g1 = [1-res_g1(:,2)./res_g1(:,1),1-res_g1(:,3)./res_g1(:,1),res_g1(:,4),find(res(:,5)==1)];
radii_relative_g2 = [1-res_g2(:,2)./res_g2(:,1),1-res_g2(:,3)./res_g2(:,1),res_g2(:,4),find(res(:,5)==2)];
radii_relative_g3 = [1-res_g3(:,2)./res_g3(:,1),1-res_g3(:,3)./res_g3(:,1),res_g3(:,4),find(res(:,5)==3)];
%% clustering group #2
clc,[idx2,C2] = kmeans(res_g2(:,1:3),2);
for ii = 1:100
    [idx2,C2] = kmeans(res_g2(:,1:3),2,"MaxIter",1000,"Start",C2);
end
colidx2 = repmat([1,0,1],length(idx2),1);
colidx2(idx2==2,:) = colidx2(idx2==2,:)*0 + [0,0.5,1];
figure(2222),scatter3(res_g2(:,1),res_g2(:,2),res_g2(:,3),30,colidx2,"filled")
% hold on
% plot(1-C2(:,2)./C2(:,1),1-C2(:,3)./C2(:,1),'kx'),drawnow;
% legend('Cluster 1','Cluster 2','Cluster Centroid')
% hold off

% gp1 = [];
% gp2 = [];
% for ii = 1:1:length(radii_relative_g2)
%     imind = radii_relative_g2(ii,4);
%     mbstack = single(tiffreadVolume(fullfile(filelist(imind).folder,filelist(imind).name)));
%     mbstackproj = rescale(max(mbstack,[],3),0,1);
%     switch idx2(ii)
%         case 1
%             gp1 = cat(2,gp1,mbstackproj);
%         case 2
%             gp2 = cat(2,gp2,mbstackproj);
%     end
%     disp("#"+ii+" "+filelist(imind).name)
% end
% figure(4445),
% subplot(211),imshow(gp1*255,hot(256)),title("[Group #1]");
% subplot(212),imshow(gp2*255,hot(256)),title("[Group #2]");
%%
clc,[idx3,C3] = kmeans(radii_relative_g3(:,1:3),2);
for ii = 1:100
    [idx3,C3] = kmeans(radii_relative_g3(:,1:3),2,"MaxIter",1000,"Start",C3);
end
colidx3 = repmat([1,0,1],length(idx3),1);
colidx3(idx3==2,:) = colidx3(idx3==2,:)*0 + [0,0.5,1];
figure(2223),scatter3(res_g3(:,1),res_g3(:,2),res_g3(:,3),30,colidx3,"filled")
% hold on
% plot(1-C3(:,2)./C3(:,1),1-C3(:,3)./C3(:,1),'kx'),drawnow;
% legend('Cluster 1','Cluster 2','Cluster Centroid')
% hold off
C3xyz = [mean(res_g3(idx3==1,1:3),1);mean(res_g3(idx3==2,1:3),1)];

% gp1 = [];
% gp2 = [];
% gp3 = [];
% for ii = 1:1:length(radii_relative_g3)
%     imind = radii_relative_g3(ii,4);
%     mbstack = single(tiffreadVolume(fullfile(filelist(imind).folder,filelist(imind).name)));
%     mbstackproj = rescale(max(mbstack,[],3),0,1);
%     switch idx3(ii)
%         case 1
%             gp1 = cat(2,gp1,mbstackproj);
%         case 2
%             gp2 = cat(2,gp2,mbstackproj);
%         case 3
%             gp3 = cat(2,gp3,mbstackproj);
%     end
%     disp("#"+ii+" "+filelist(imind).name)
% end
% figure(4445),
% subplot(311),imshow(gp1*255,hot(256)),title("[Group #1]");
% subplot(312),imshow(gp2*255,hot(256)),title("[Group #2]");
% subplot(313),imshow(gp3*255,hot(256)),title("[Group #3]");

%% --------------------------------------------------------------
% res_g31 = res_g3(idx3==1,:);
% radii_relative_g31 = radii_relative_g3(idx3==1,:);
% clc,[idx31,C31] = kmeans(res_g31(:,1:3),2);
% for ii = 1:100
%     [idx31,C31] = kmeans(res_g31(:,1:3),2,"MaxIter",1000,"Start",C31);
% end
% colidx31 = repmat([1,0,1],length(idx31),1);
% colidx31(idx31==2,:) = colidx31(idx31==2,:)*0 + [0,0.5,1];
% figure(2222),scatter3(res_g31(:,1),res_g31(:,2),res_g31(:,3),30,colidx31,"filled")
% % hold on
% % plot(1-C3(:,2)./C3(:,1),1-C3(:,3)./C3(:,1),'kx'),drawnow;
% % legend('Cluster 1','Cluster 2','Cluster Centroid')
% % hold off
% 
% gp1 = [];
% gp2 = [];
% for ii = 1:1:length(radii_relative_g31)
%     imind = radii_relative_g31(ii,4);
%     mbstack = single(tiffreadVolume(fullfile(filelist(imind).folder,filelist(imind).name)));
%     mbstackproj = rescale(max(mbstack,[],3),0,1);
%     switch idx3(ii)
%         case 1
%             gp1 = cat(2,gp1,mbstackproj);
%         case 2
%             gp2 = cat(2,gp2,mbstackproj);
%         case 3
%             gp3 = cat(2,gp3,mbstackproj);
%     end
%     disp("#"+ii+" "+filelist(imind).name)
% end
% figure(4445),
% subplot(311),imshow(gp1*255,hot(256)),title("[Group #1]");
% subplot(312),imshow(gp2*255,hot(256)),title("[Group #2]");
% subplot(313),imshow(gp3*255,hot(256)),title("[Group #3]");
%% ========================================================
clc,res_all = res;
res_all(radii_relative_g2(idx2==1,4),5) = 2;
res_all(radii_relative_g2(idx2==2,4),5) = 3;
res_all(radii_relative_g3(idx3==1,4),5) = 4;
res_all(radii_relative_g3(idx3==2,4),5) = 5;
colidx = repmat([1,0,1],length(radii_relative),1); % connected two cells
colidx(res_all(:,5)==2,:) = colidx(res_all(:,5)==2,:)*0 + [1,0.5,0.2]; 
colidx(res_all(:,5)==3,:) = colidx(res_all(:,5)==3,:)*0 + [0,0.7,0.1];
colidx(res_all(:,5)==4,:) = colidx(res_all(:,5)==4,:)*0 + [0.1,0.8,0.8];
colidx(res_all(:,5)==5,:) = colidx(res_all(:,5)==5,:)*0 + [0.5,0.5,1];

fig = figure(4444);
set(fig, 'Position',  [600, 200, 1080, 600])
ax = axes('Parent',fig); %view(ax,[-40 35]);hold(ax,'off');
scattersize = 64;
scatter3(res_all(res_all(:,5)==1,1),res_all(res_all(:,5)==1,2),res_all(res_all(:,5)==1,3),...
    scattersize,colidx(res_all(:,5)==1,:),"filled"),hold on;
scatter3(res_all(res_all(:,5)==2,1),res_all(res_all(:,5)==2,2),res_all(res_all(:,5)==2,3),...
    scattersize,colidx(res_all(:,5)==2,:),"filled"),hold on;
scatter3(res_all(res_all(:,5)==3,1),res_all(res_all(:,5)==3,2),res_all(res_all(:,5)==3,3),...
    scattersize,colidx(res_all(:,5)==3,:),"filled"),hold on;
scatter3(res_all(res_all(:,5)==4,1),res_all(res_all(:,5)==4,2),res_all(res_all(:,5)==4,3),...
    scattersize,colidx(res_all(:,5)==4,:),"filled"),hold on;
scatter3(res_all(res_all(:,5)==5,1),res_all(res_all(:,5)==5,2),res_all(res_all(:,5)==5,3),...
    scattersize,colidx(res_all(:,5)==5,:),"filled"),hold on;
xlabel("R_a (μm)"),ylabel("R_b (μm)"),zlabel("R_c (μm)");
xlim([30 95]),xticks(30:15:90)
ylim([20 40]),yticks(20:5:40)
zlim([5 40]),zticks(5:7:40)
set(ax,"fontsize",28)
set(ax,'Color','none','FontSize',28,'GridAlpha',1,'GridColor',[0.8 0.8 0.8])
% legend("Adhered cells","Small spheres","Big spheres","Sickle/Broken cells","Ovalcytes/Helmet cells")
%%
xsize = 1080;
ysize = 460;
% plot scatter on xy plane
fig = figure(51);
set(fig, 'Position',  [600, 200, xsize, ysize])
ax = axes('Parent',fig); %view(ax,[-40 35]);hold(ax,'off');
scattersize = 64;
scatter(res_all(res_all(:,5)==1,1),res_all(res_all(:,5)==1,2),...
    scattersize,colidx(res_all(:,5)==1,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==2,1),res_all(res_all(:,5)==2,2),...
    scattersize,colidx(res_all(:,5)==2,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==3,1),res_all(res_all(:,5)==3,2),...
    scattersize,colidx(res_all(:,5)==3,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==4,1),res_all(res_all(:,5)==4,2),...
    scattersize,colidx(res_all(:,5)==4,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==5,1),res_all(res_all(:,5)==5,2),...
    scattersize,colidx(res_all(:,5)==5,:),"filled"),hold on;
xlabel("R_a (μm)"),ylabel("R_b (μm)")
xlim([30 95]),xticks(30:15:90)
ylim([20 40]),yticks(20:5:40)
set(ax,"xGrid",true,"yGrid",true)
set(ax,"Box",true)
set(ax,"fontsize",28)
set(ax,'Color','none','FontSize',28,'GridAlpha',1,'GridColor',[0.8 0.8 0.8])
%
% plot scatter on xz plane
fig = figure(52);
set(fig, 'Position',  [600, 200, xsize, ysize])
ax = axes('Parent',fig); %view(ax,[-40 35]);hold(ax,'off');
scattersize = 64;
scatter(res_all(res_all(:,5)==1,1),res_all(res_all(:,5)==1,3),...
    scattersize,colidx(res_all(:,5)==1,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==2,1),res_all(res_all(:,5)==2,3),...
    scattersize,colidx(res_all(:,5)==2,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==3,1),res_all(res_all(:,5)==3,3),...
    scattersize,colidx(res_all(:,5)==3,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==4,1),res_all(res_all(:,5)==4,3),...
    scattersize,colidx(res_all(:,5)==4,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==5,1),res_all(res_all(:,5)==5,3),...
    scattersize,colidx(res_all(:,5)==5,:),"filled"),hold on;
xlabel("R_a (μm)"),ylabel("R_c (μm)");
xlim([30 95]),xticks(30:15:90)
ylim([5 40]),yticks(5:7:40)
set(ax,"xGrid",true,"yGrid",true)
set(ax,"Box",true)
set(ax,"fontsize",28)
set(ax,'Color','none','FontSize',28,'GridAlpha',1,'GridColor',[0.8 0.8 0.8])
%
% plot scatter on yz plane
fig = figure(53);
set(fig, 'Position',  [600, 200, xsize, ysize])
ax = axes('Parent',fig); %view(ax,[-40 35]);hold(ax,'off');
scattersize = 64;
scatter(res_all(res_all(:,5)==1,2),res_all(res_all(:,5)==1,3),...
    scattersize,colidx(res_all(:,5)==1,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==2,2),res_all(res_all(:,5)==2,3),...
    scattersize,colidx(res_all(:,5)==2,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==3,2),res_all(res_all(:,5)==3,3),...
    scattersize,colidx(res_all(:,5)==3,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==4,2),res_all(res_all(:,5)==4,3),...
    scattersize,colidx(res_all(:,5)==4,:),"filled"),hold on;
scatter(res_all(res_all(:,5)==5,2),res_all(res_all(:,5)==5,3),...
    scattersize,colidx(res_all(:,5)==5,:),"filled"),hold on;
xlabel("R_b (μm)"),ylabel("R_c (μm)");
set(ax,"Box",true)
xlim([20 40]),xticks(20:5:40)
ylim([5 40]),yticks(5:7:40)
set(ax,"xGrid",true,"yGrid",true)
set(ax,"fontsize",28)
set(ax,'Color','none','FontSize',28,'GridAlpha',1,'GridColor',[0.8 0.8 0.8])