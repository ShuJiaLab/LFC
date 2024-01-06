clear,clc

rawtif_folder = "Z:\Xuanwen\FLFMuf\ExpData\XW20220804_VLS\V6\";
if ~exist(rawtif_folder + "rawtif_expressed_a",'dir')
    mkdir(rawtif_folder + "rawtif_expressed_a");
end
if ~exist(rawtif_folder + "rawtif_expressed_t",'dir')
    mkdir(rawtif_folder + "rawtif_expressed_t");
end
if ~exist(rawtif_folder + "rawtif_normal_a",'dir')
    mkdir(rawtif_folder + "rawtif_normal_a");
end
%% decide which frame is membrane / cytoplasm
% Rename / switch the folder names
% This should be done manually
% V5 has been verified for a/t.

% folder1 = "rawtif_1_a";
% folder2 = "rawtif_1_t";
% movefile(folder1,"rawtif_temp")
% movefile(folder2,folder1)
% movefile("rawtif_temp",folder2)
%% Get filenames and fileindeces
filenames = dir(rawtif_folder+ "rawtif" +"\*.tif");
lfprefix = "V6";
filelist = {filenames(:).name}';
fileinds = split(filelist(:),["(",")","_",".tif"]);
filelist = [fileinds(:,1),num2cell(str2double(fileinds(:,[2,3,5]))),filelist];
h = waitbar(0,'Getting file page numbers...');
for ii = 1:1:length(filelist)
    filelist{ii,4} = length(imfinfo(rawtif_folder + "rawtif\" + filelist{ii,5}));
    waitbar(ii/length(filelist),h,'Getting file page numbers...');
end
filelist = sortrows(filelist);
close(h)
filelist(string(filelist(:,1))~=lfprefix,:)=[];
%% ========================================================================
% For any non-verified a/t groups, the default sequence is odd for "a" and
% even for "t".
oddisa = [0,1,1,1,1,           1,1,1,1,1];
threslow = 200;
threshigh = 65500;
baseline = 90;
% =========================================================================
%% Get filenames and fileindeces (OBSOLETE)
% groupname = 'rawtif_1_a';
% expgroup = groupname(1:end-1) + "t";
% filenames = dir(rawtif_folder+groupname+"\*.tif");
% filelist = {filenames(:).name}';
% fileinds = split(filelist(:),["(",")","_",".tif"]);
% filelist = [fileinds(:,1),num2cell(str2double(fileinds(:,[2,3,5]))),filelist];
% filelist = sortrows(filelist);
%% Select 100
clc,
try
    load(rawtif_folder+"readpt.mat","readpt");
catch
    readpt = [2,1,1,1];
    save(rawtif_folder+"readpt.mat","readpt");
end
% ------------------------------------------
groupind = readpt(1);
subgroupind = readpt(2);
pageind = readpt(3);
totalpageind = readpt(4);
for ii = 2%groupind:1:filelist{end,2}
    for jj = subgroupind:1:max(cell2mat(filelist(cell2mat(filelist(:,2))==ii,3)))
        filename = filelist{logical((cell2mat(filelist(:,2))==ii) .* (cell2mat(filelist(:,3))==jj)),5};
        for kk = pageind:1:filelist{string(filelist(:,5))==string(filename),4}
            if mod(totalpageind,2) == oddisa(ii)
                try
                    im = imread(rawtif_folder+ "rawtif\" + filename,'Index',kk);
                catch
                    im = ones(900,1024,3);
                    im(1,1)=0;
                end
                if max(im(:))>=threslow && max(im(:))<=threshigh
                    try
                        im1 = imread(rawtif_folder+ "rawtif\" + filename,'Index',kk-1);
                    catch
                        im1 = ones(900,1024,3);
                        im1(1,1)=0;
                    end
                    try
                        im2 = imread(rawtif_folder+"rawtif\" + filename,'Index',kk+1);
                    catch
                        im2 = ones(900,1024,3);
                        im2(1,1)=0;
                    end
                    imbatch = cat(2,...
                        insertText(rescale(single(im1-baseline),0,1),[700,800],...
                        filelist{ii,1} + "_" + ii + "(" + jj + ")" + "_" + (kk-1),...
                        'Font','Calibri Bold','FontSize',48,'TextColor','green','BoxOpacity',0),...
                        insertText(rescale(single(im-baseline),0,1),[700,800],...
                        filelist{ii,1} + "_" + ii + "(" + jj + ")" + "_" + kk,...
                        'Font','Calibri Bold','FontSize',48,'TextColor','red','BoxOpacity',0),...
                        insertText(rescale(single(im2-baseline),0,1),[700,800],...
                        filelist{ii,1} + "_" + ii + "(" + jj + ")" + "_" + (kk+1),...
                        'Font','Calibri Bold','FontSize',48,'TextColor','green','BoxOpacity',0) );
                    savefilename = rawtif_folder + "rawtif_normal_a\" + filename(1:end-4)+"_"+kk+".tif";
                    savefilename1 = rawtif_folder + "rawtif_expressed_a\" + filename(1:end-4)+"_"+kk+".tif";
                    savefilename2 = rawtif_folder + "rawtif_expressed_t\" + filename(1:end-4)+"_"+(kk-1)+".tif";
                    savefilename3 = rawtif_folder + "rawtif_expressed_t\" + filename(1:end-4)+"_"+(kk+1)+".tif";
                    colsize = 3;
                    [row,col,~] = size(imbatch);
                    col = round(col/colsize);
                    try
                        batchlist = [savefilename,savefilename1,savefilename2,savefilename3];
                        imbatch = updatebatch(imbatch,batchlist,row,col);
                        imgtitle = filelist{ii,1}+" "+ii+" ("+jj+") "+...
                            (kk-1)+'-'+ (kk+1) + " | "+ round(ii/length(filelist)*100,2)+"%";
                        imorder = batch2order(imbatch,row,col,colsize,imgtitle);
                    catch
                        readpt = [ii,jj,kk,totalpageind];
                        save(rawtif_folder+"readpt.mat","readpt");
                        return
                    end
                    if length(imorder)==1
                        imwrite(im,savefilename);
                        disp(filename(1:end-4)+"_"+ kk + ' normal')
                    elseif length(imorder)==2
                        imwrite(im,savefilename1); % 1st click is membrance
                        if imorder(2)==1           % 2nd click is cytoplasma
                            imwrite(im1,savefilename2);
                            disp(filename(1:end-4)+"_"+ kk + ' expressed - axlexa647 | '+ ...
                                filename(1:end-4)+"_"+ (kk-1) +' expressed - tdtomato');
                        elseif imorder(2)==3
                            imwrite(im2,savefilename3);
                            disp(filename(1:end-4)+"_"+ kk + ' expressed - axlexa647 | '+ ...
                                filename(1:end-4)+"_"+ (kk+1) +' expressed - tdtomato');
                        end
                    else
                        disp(filename(1:end-4)+"_"+ kk + "[No cell selected!]")
                    end
                else
                    disp(filelist{ii,1}+" "+ii+" ("+jj+") "+...
                        (kk-1)+'-'+ (kk+1) + " | "+ round(ii/length(filelist)*100,2)+"%" + " [TOO DIM]");
                end
            end
            readpt = [ii,jj,kk,totalpageind];
            save(rawtif_folder+"readpt.mat","readpt");
            totalpageind = totalpageind + 1;
        end
        pageind = 1; % reset pageind
    end
    subgroupind = 1; % reset subgroupind
    totalpageind = 1; % reset totalpageind
end

%% ====================== Utilities ========================
% Here contains all the subfunctions that are used in the main function.
% ==========================================================

function imorder = batch2order(imbatch,row,col,colsize,imgtitle)
figure(99),imshow(imbatch),title(imgtitle);
[posx,posy] = ginput(2);
if length(posx)==2
    imorder = [floor(posy(1)/row)*colsize+ceil(posx(1)/col),...
               floor(posy(2)/row)*colsize+ceil(posx(2)/col)];
elseif length(posx)==1
    imorder = floor(posy(1)/row)*colsize+ceil(posx(1)/col);
else
    imorder = [];
end
close(figure(99))
end

function imbatch = updatebatch(imbatch,filelist,filerow,filecol)
if exist(filelist(1),"file") || ...
   exist(filelist(2),"file")
    imbatch(1:filerow,filecol+1:2*filecol,:) = 1;
end
if exist(filelist(3),"file")
    imbatch(1:filerow,1:filecol,:) = 1;
end
if exist(filelist(4),"file")
    imbatch(1:filerow,2*filecol+1:3*filecol,:) = 1;
end
end