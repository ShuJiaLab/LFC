clear,clc

params = ["XW20220609","XW20220616","XW20220618","XW20220619","XW20220620";
    "JK_minu","JK_minu30_","JK_minu60_","JK_minu120_","JK_minu300_"];

groupind = 5;
groupname = params(1,groupind) ;
prefix = params(2,groupind);
nufolder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220616_JurkatSTS\"+groupname+"\PSFFLFint_20220630_Blue_Gly_10um-z_multi_nucleus_acsn_ccut";
nufolderlabeled = nufolder + "_labeled\";
filenames = dir(nufolderlabeled+"\*.tif");
filelist = {filenames(:).name}';
fileinds = split(filelist(:),[prefix,"(",")","_","iter50.tif"]);
filelist = [fileinds(:,1),num2cell(str2double(fileinds(:,[2,3,5]))),filelist];
filelist = sortrows(filelist);
filelist(:,1)={char(prefix)};
%%
NusomeDist = [];
for ii =1:1:length(filelist) 
    if ~exist(nufolderlabeled+"\"+filelist{ii,5},"file")
        return
    else
        imstacklabeled = tiffreadVolume(nufolderlabeled+"\"+filelist{ii,5});
        numofnu = max(imstacklabeled(:));
        vols = zeros(numofnu,2);
        for jj = 1:1:numofnu
            vols(jj,:) = [jj,sum(double(imstacklabeled==jj),"all")];
        end
        vols(:,2) = double(vols(:,2))*0.065*0.065*0.1;
        vols(vols(:,2)<pi/6*0.9^3,:)=[];
        ind = [];
        for kk = 1:1:size(vols,1)
            ind = [ind;find(imstacklabeled == vols(kk,1))];
        end
        [rows,cols,pages] = ind2sub(size(imstacklabeled),ind);
        cent = mean([rows,cols,pages],1);
        dist = [];
        for kk = 1:1:size(vols,1)
            [rows,cols,pages] = ind2sub(size(imstacklabeled),find(imstacklabeled == vols(kk,1)));
            centi = mean([rows,cols,pages],1);
            dist = [dist,sqrt(sum((centi-cent).^2,"all"))];
        end
        NusomeDist = cat(1,NusomeDist,{filelist{ii,5},numofnu,dist,mean(dist),std(dist)*(length(dist)-1)});
        disp(ii+" ["+filelist{ii,5}+"] "+numofnu+": ["+num2str(dist)+"] * Saved! *")
    end
end
%%
clc;
dist = [];
for ii =1:1:length(filelist) 
    if ~exist(nufolderlabeled+"\"+filelist{ii,5},"file")
        return
    else
        imstacklabeled = tiffreadVolume(nufolderlabeled+"\"+filelist{ii,5});
        numofnu = max(imstacklabeled(:));
        vols = zeros(numofnu,2);
        for jj = 1:1:numofnu
            vols(jj,:) = [jj,sum(double(imstacklabeled==jj),"all")];
        end
        vols(:,2) = double(vols(:,2))*0.065*0.065*0.1;
        vols(vols(:,2)<pi/6*0.9^3,:)=[];
        ind = [];
        for kk = 1:1:size(vols,1)
            ind = [ind;find(imstacklabeled == vols(kk,1))];
        end
        [rows,cols,pages] = ind2sub(size(imstacklabeled),ind);
        cent = mean([rows,cols,pages],1);
        
        for kk = 1:1:size(vols,1)
            [rows,cols,pages] = ind2sub(size(imstacklabeled),find(imstacklabeled == vols(kk,1)));
            centi = mean([rows,cols,pages],1);
            dist = [dist,sqrt(sum((centi-cent).^2,"all"))];
        end
        disp(ii+" ["+filelist{ii,5}+"] "+numofnu+": ["+num2str(length(dist))+"] * Saved! *")
    end
end
dist = dist';