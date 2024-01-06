clear,clc

foldername = "./PSFFLFint_20220515_Blue_Gly_10um-rawtif_selected_focused_offsetsub_acsn_ccut/";
filelist = dir(foldername+"/*.tif");
%%
clc,res = zeros(100,4);
%%
for ii = 100 %:length(filelist)
mbstack = single(tiffreadVolume(fullfile(filelist(ii).folder,filelist(ii).name)));
mbstack = imresize3(mbstack,round([size(mbstack,1),size(mbstack,2),size(mbstack,3)*100/65]));
mbstack = mbstack./max(mbstack,[],"all");
% mbstack(127-30:127+30,190-30:190+10,:) = min(mbstack(:));
% mbstack(100-30:100+30,33-30:33+30,:) = min(mbstack(:));
fudgeFactor = 0.4;
imBWsstack = edge3(double(mbstack>0.4*fudgeFactor).*double(mbstack),'sobel',0.5*fudgeFactor);
imBWsstackfilled = imclose(imBWsstack,strel('disk',30));
figure(1),sliceViewer([mbstack,single(imBWsstackfilled)]);
[row,col,depth] = ind2sub(size(imBWsstackfilled),find(imBWsstackfilled>0));
xyzPoints = [row,col,depth];
[ center, radii, ~, ~, chi2 ] = ellipsoid_fit( xyzPoints, '' );
V = sum(imBWsstackfilled,"all")*0.065^3;
disp("[#"+filelist(ii).name+"] radii: "+radii(1)+", "+radii(2)+", "+radii(3)+"; V: "+V)
res(ii,:) = [radii' ,V];
% figure(1),sliceViewer([mbstack,single(imBWsstackfilled)]);
end
save("res.mat","res");
%%
clc,
radii_relative = [1-res(:,2)./res(:,1),1-res(:,3)./res(:,1),res(:,4)];
% radii_relative(radii_relative(:,1)<0,:)=[];
e1cutoff = 1-1/1.35;
e2cutoff1 = 1-1/1.42;
e2cutoff2 = 1-1/2.55;
Vcutoff1 = 60;
Vcutoff2 = 110;
divideN = 800;
groupN = 16;
[~,Xedges,Yedges] = histcounts2(radii_relative(:,1),radii_relative(:,2),divideN);
[N,~,~] = histcounts2(radii_relative(:,1),radii_relative(:,2),groupN);
% [N,~,~] = histcounts2(log(radii_relative(:,1)),log(radii_relative(:,2)),groupN);
N = imresize(N,round(divideN/groupN));
figure(1),imagesc(flipud(N'));

Cvals = [];
Nvals = [];
Svals = [];
figure(2),
for ii = 1:1:length(radii_relative)
    if radii_relative(ii,3)<Vcutoff1
        Svalue = "o";
    elseif radii_relative(ii,3)>Vcutoff2
        Svalue = "+";
    else
        Svalue = "filled";
    end

    Nvalue = getNvalue(radii_relative(ii,1),radii_relative(ii,2),Xedges,Yedges,N);
%     Cvalue = getCvalue(cool(256),Nvalue,min(N(:)),max(N(:)));
    if radii_relative(ii,1)<(e1cutoff)
        if radii_relative(ii,2)<(e2cutoff1)
            Cvalue = [1,0.5,0.5];
        else
            Cvalue = [1,0.7,0.1];
        end
    else
        if radii_relative(ii,2)<(e2cutoff2)
            Cvalue = [0,0.8,0.8];
        else
            Cvalue = [0.2,0.5,1];
        end
    end
    if ~ismember(ii,[44,43,77,26,90,54])
        scatter3(radii_relative(ii,1),radii_relative(ii,2),radii_relative(ii,3),...
                                                    56,Cvalue,Svalue,"LineWidth",2);hold on;
    else
        scatter3(radii_relative(ii,1),radii_relative(ii,2),radii_relative(ii,3),...
                                                    56,[0,0,0],"*","LineWidth",2);hold on;
    end
    drawnow;
    Cvals = cat(1,Cvals,Cvalue);
    Nvals = cat(1,Nvals,Nvalue);
    Svals = cat(1,Svals,Svalue);
end
hold off;
% scatter3(radii_relative(:,1),radii_relative(:,2),radii_relative(:,3),64,Cvals,"o");
% set(gca,'xscale','log')
% set(gca,'yscale','log')
%%
clc,
radii_relative_temp = radii_relative;
figure(100)
for ii = 1:1:100
    row = find(radii_relative_temp(:,1)==min(radii_relative_temp(:,1)));
    mbstack = single(tiffreadVolume(fullfile(filelist(row).folder,filelist(row).name)));
    mbstackproj = max(mbstack,[],3);
    mbstackproj = mbstackproj/max(mbstackproj(:));
    subplot(10,10,ii),imshow(mbstackproj*255,hot(256)),
    title("["+round(radii_relative(row,1),2)+","+round(radii_relative(row,2),2)+","+round(radii_relative(row,3),2)+"]"),
    drawnow;
    radii_relative_temp(row,1)=100;
    disp("#"+ii+" "+row+" "+filelist(row).name)
end
%%
clc,
radii_relative_temp = radii_relative;
figure(400)
e1cutoff = 1-1/1.35;
e2cutoff1 = 1-1/1.42;
e2cutoff2 = 1-1/2.55;
Vcutoff1 = 20;
Vcutoff2 = 40;
e1se2s = [];
e1se2l = [];
e1le2s = [];
e1le2l = [];
for ii = 1:1:100
    mbstack = single(tiffreadVolume(fullfile(filelist(ii).folder,filelist(ii).name)));
    mbstackproj = rescale(max(mbstack,[],3),0,1);
    if radii_relative_temp(ii,1)<=e1cutoff
        if radii_relative_temp(ii,2)<=e2cutoff1
            e1se2s = cat(2,e1se2s,mbstackproj);
        else
            e1se2l = cat(2,e1se2l,mbstackproj);
        end
    else
        if radii_relative_temp(ii,2)<=e2cutoff2
            e1le2s = cat(2,e1le2s,mbstackproj);
        else
            e1le2l = cat(2,e1le2l,mbstackproj);
        end
    end
    subplot(411),imshow(e1se2s*255,hot(256)),title("[e1 small e2 small]"),drawnow;
    subplot(412),imshow(e1se2l*255,hot(256)),title("[e1 small e2 large]"),drawnow;
    subplot(413),imshow(e1le2s*255,hot(256)),title("[e1 large e2 small]"),drawnow;
    subplot(414),imshow(e1le2l*255,hot(256)),title("[e1 large e2 large]"),drawnow;
    disp("#"+ii+" "+filelist(row).name)
end

%%
res = res./max(res,[],1);
[coeff,score,latent]=pca(res); 
latentnorm = cumsum(latent)/sum(latent);
biplot(coeff(:,1:3),'Scores',score(:,1:3),'VarLabels',{'Ra' 'Rb' 'Rc' 'Chi2'});
df=score(:,1:3);
y=pdist(df,'euclidean');
z=linkage(y,'average');
T=cluster(z,'maxclust',3);
for i=1:3
     tm=find(T==i); 
     na='';
      for j=1:length(tm)
          na=strcat(na,',',num2str(tm(j)));
      end
       fprintf('第%d类的有%s\n',i,na);
end

%% K-means clustering
foldername = "./PSFFLFint_20220515_Blue_Gly_10um-rawtif_selected_focused_offsetsub_acsn_ccut/";
filelist = dir(foldername+"/*.tif");
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

e1se2s = [];
e1se2l = [];
e1le2s = [];
e1le2l = [];
for ii = 1:1:100
    mbstack = single(tiffreadVolume(fullfile(filelist(ii).folder,filelist(ii).name)));
    mbstackproj = rescale(max(mbstack,[],3),0,1);
    switch idx(ii)
        case 1
            e1se2s = cat(2,e1se2s,mbstackproj);
        case 2
            e1se2l = cat(2,e1se2l,mbstackproj);
        case 3
            e1le2s = cat(2,e1le2s,mbstackproj);
%         case 4
%             e1le2l = cat(2,e1le2l,mbstackproj);
    end
    disp("#"+ii+" "+filelist(ii).name)
end
figure(4445),
subplot(411),imshow(e1se2s*255,hot(256)),title("[Group #1]");
subplot(412),imshow(e1se2l*255,hot(256)),title("[Group #2]");
subplot(413),imshow(e1le2s*255,hot(256)),title("[Group #3]");
subplot(414),imshow(e1le2l*255,hot(256)),title("[Group #4]");
%% utilities

function Nvalue = getNvalue(beadvol,beadint,Xedges,Yedges,N)
Xpos = length(find(beadvol>=Xedges));
Ypos = length(find(beadint>=Yedges));
% disp([Xpos,Ypos]);
Nvalue = N(Xpos,Ypos);
end

function Cvalue = getCvalue(cmp,Nvalue,Nmin,Nmax)
Cind = round(255/(Nmax-Nmin)*(Nvalue-Nmin)+1);
Cvalue = cmp(Cind,:);
end