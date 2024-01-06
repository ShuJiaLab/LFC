clear,clc

thres = 350;
addpath('..\utilities\')
scrprefix = '.\rawtif\pemi';

frameRes = [];
H = multiwaitbar(4,[0 0 0 0],{'Reading Group Pages...','Reading Sub-Groups...',...
                              'Reading Exp Groups...','Frames collected: 0'});
frameCt = 0;
ii_start = 1;
ii_max = 9;
jj_max = 26;
cell_max = 1e5;
for ii = ii_start:1:ii_max
    for jj = 1:1:jj_max
        scrname = [scrprefix,num2str(ii,'%05d'),'(',num2str(jj),').tif'];
%         disp(scrname)
        % ===================
        kk = 1;
        while 1
            try
                im = imread(scrname,'index',kk);
                immax = max(double(im(:)));
                if immax>thres
                    frameRes = cat(1,frameRes,...
                        [ii,jj,kk,immax,sum(im>thres,'all')]);
                    frameCt = frameCt+1;
                end
                multiwaitbar(4,...
                    [((jj-1)*2329+kk)/60000,jj/jj_max,(ii-1)/(ii_max-ii_start+1),frameCt/cell_max],...
                    {['Reading Page #',num2str((jj-1)*2329+kk),'/60000'],...
                    ['Reading Sub-Group #',num2str(jj),'/',num2str(jj_max)],...
                    ['Reading Exp Group #',num2str(ii),'/',num2str(ii_start)],...
                    ['Frames collected: ',num2str(frameCt),'/',num2str(cell_max)]},...
                    H);

                kk = kk + 1;
            catch ME
                if kk==1
                    rethrow(ME)
                else
                    break
                end
            end
        end
        % ===================
    end
end
delete(H.figure)
save("frameRes.mat","frameRes")
%% Brief selection (bright area size)

load("frameRes.mat","frameRes")
frameRes1 = frameRes;
frameRes1(frameRes1(:,5)<100,:)=[];
%% Fine selection (continuous frames)
frameRes2 = [];
H = multiwaitbar(1,[0],{'ind: 0'});
for ind = 1:1:(size(frameRes1,1)-1)
    if ind == 1
        ctn = (frameRes1(ind+1,3) == (frameRes1(ind,3)+1));
    else
        ctn = (frameRes1(ind+1,3) == (frameRes1(ind,3)+1) | frameRes1(ind-1,3) == (frameRes1(ind,3)-1));
    end
    if ctn
        frameRes2 = cat(1,frameRes2,frameRes1(ind,:));
    end
    multiwaitbar(1,[ind/(size(frameRes1,1)-1)],{['ind: ',num2str(ind),'/',num2str(size(frameRes1,1)-1)]},H);
end
delete(H.figure)
save("frameRes2.mat","frameRes2")
%% Cell counting
clc,
load("frameRes2.mat","frameRes2")
cellcount = [];
H = multiwaitbar(1,[0],{'cellcount: []'});
for ind = 1:1:(size(frameRes2,1)-1)
    if ind == 1
        ctn = (frameRes2(ind+1,3) == frameRes2(ind,3));
        cellcount = [frameRes2(ind,1),frameRes2(ind,2),frameRes2(ind,3),frameRes2(ind,3),0];
        multiwaitbar(1,[ind/(size(frameRes2,1)-1)],{['#',num2str(ind),': ',num2str(cellcount(end,:))]},H);
    else
        ctn = (frameRes2(ind+1,3) == (frameRes2(ind,3)+1));
        if ~ctn
            cellcount(end,4) = frameRes2(ind,3);
            cellcount(end,5) = cellcount(end,4) - cellcount(end,3) + 1;
            cellcount = cat(1,cellcount,...
                [frameRes2(ind+1,1),frameRes2(ind+1,2),frameRes2(ind+1,3),frameRes2(ind+1,3),0]);
            multiwaitbar(1,[ind/(size(frameRes2,1)-1)],{['#',num2str(ind),': ',num2str(cellcount(end-1,:))]},H);
        end
    end
end
% cellcount(end,4) = frameRes2(ind,3);
delete(H.figure)
% cellcount(end,5) = cellcount(end,4) - cellcount(end,3) + 1;
save("cellcount.mat","cellcount")
%%
clc
addpath('..\utilities\')
load('frameRes2.mat','frameRes2')
thres_2 = 350;
if ~exist('rawtif_selected','dir')
    mkdir('.\rawtif_selected')
end
H = multiwaitbar(1,[0],{'Selecting cells...'});
for cellind = 1:1:size(frameRes2,1)
    im1 = imread(['rawtif\pemi',num2str(frameRes2(cellind,1),'%05d'),'(',...
        num2str(frameRes2(cellind,2)),').tif'],"index",frameRes2(cellind,3));
    im1 = [repmat(im1(1,:),62,1);im1;repmat(im1(end,:),62,1)];
    if (max(im1(:))>thres_2)% && max(im2(:))>thres_2)
        imwrite(im1,['rawtif_selected\pemi',num2str(frameRes2(cellind,1),'%05d'),...
            '(',num2str(frameRes2(cellind,2)),')_',num2str(frameRes2(cellind,3)),'.tif']);
    multiwaitbar(1,[cellind/size(frameRes2,1)],{['Selected (',num2str(round(cellind/size(frameRes2,1)*100)),'%)']},H)
    end
end
delete(H.figure)
%%
clc
addpath('..\utilities\')
load('cellcount.mat','cellcount')
thres_2 = 350;
if ~exist('rawtif_selected_pe','dir')
    mkdir('rawtif_selected_pe')
    mkdir('rawtif_selected_mi')
end
H = multiwaitbar(1,[0],{'Selecting cells...'});
cellselected = [];
for cellind = 918:1:size(cellcount,1)
    cellcountext = [repmat(cellcount(cellind,1:2),cellcount(cellind,5),1),(cellcount(cellind,3):1:cellcount(cellind,4))'];
    cellfilelist = [];
    for ii = 1:1:cellcount(cellind,5)
        cellfilename = string(['rawtif_selected\pemi',num2str(cellcountext(ii,1),'%05d'),'(',...
        num2str(cellcountext(ii,2)),')_',num2str(cellcountext(ii,3)),'.tif']);
        cellfilelist = cat(1,cellfilelist,cellfilename);
    end
    cells = cellctext(cellcountext,cellfilelist,2);
    if length(cells)>2
        prompt = {'Cell 1:','Cell 2:'};
        dlgtitle = ['[',num2str(cellcount(cellind,1:4)),']'];
        dims = [1 50];
        definput = {num2str(cellcount(cellind,3)),num2str(cellcount(cellind,4))};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer)
            answer2 = questdlg('Would you like to skip this cell?', ...
                'Skip Cell', 'Skip','Cancel','Skip');
            switch answer2
                case 'Skip'
                    continue
                case 'Cancel'
                    prompt = {'Cell 1:','Cell 2:'};
                    dlgtitle = ['[',num2str(cellcount(cellind,1:4)),']'];
                    dims = [1 50];
                    definput = {num2str(cellcount(cellind,3)),num2str(cellcount(cellind,4))};
                    answer = inputdlg(prompt,dlgtitle,dims,definput);
            end
        end
        cells = [str2double(answer{1}),str2double(answer{2})];
    end
    cell_1 = cells(1);
    cell_2 = cells(2);
    im1 = imread(['rawtif_selected\pemi',num2str(cellcount(cellind,1),'%05d'),'(',...
        num2str(cellcount(cellind,2)),')_',num2str(cell_1),'.tif']);
    im2 = imread(['rawtif_selected\pemi',num2str(cellcount(cellind,1),'%05d'),'(',...
        num2str(cellcount(cellind,2)),')_',num2str(cell_2),'.tif']);

    if (max(im1(:))>thres_2 && max(im2(:))>thres_2)
        
        im1_cp = double(im1)/max(double(im1(:)));
        im2_cp = double(im2)/max(double(im2(:)));
        t1 = thres_2/double(max(im1(:)));
        t2 = thres_2/double(max(im2(:)));
        t = max([t1,t2]);

        se = strel('disk',100);
        im1_bw = imbinarize(im1_cp,t);
        im1_closed = imclose(im1_bw,se);
        sum1 = sum(im1_closed-im1_bw,'all');
%             figure(1),imshowpair(im1_bw,im1_closed,'montage');
        im2_bw = imbinarize(im2_cp,t);
        im2_closed = imclose(im2_bw,se);
        sum2 = sum(im2_closed-im2_bw,'all');
%             figure(2),imshowpair(im2_bw,im2_closed,'montage');
        if sum1>sum2
            cellselected = cat(1,cellselected,[cellcount(cellind,1),cellcount(cellind,2),cell_2,cell_1]);
            imwrite(im1,['rawtif_selected_mi\pemi',num2str(cellcount(cellind,1),'%05d'),'(',num2str(cellcount(cellind,2)),')_',num2str(cell_1),'.tif']);
            imwrite(im2,['rawtif_selected_pe\pemi',num2str(cellcount(cellind,1),'%05d'),'(',num2str(cellcount(cellind,2)),')_',num2str(cell_2),'.tif']);
            multiwaitbar(1,[cellind/size(cellcount,1)],{['Selected (',num2str(round(cellind/size(cellcount,1)*100)),'%): ',...
                num2str(cell_2),' (Perox), ',num2str(cell_1),' (Mitok).']},H)
        else
            cellselected = cat(1,cellselected,[cellcount(cellind,1),cellcount(cellind,2),cell_1,cell_2]);
            imwrite(im1,['rawtif_selected_pe\pemi',num2str(cellcount(cellind,1),'%05d'),'(',num2str(cellcount(cellind,2)),')_',num2str(cell_1),'.tif']);
            imwrite(im2,['rawtif_selected_mi\pemi',num2str(cellcount(cellind,1),'%05d'),'(',num2str(cellcount(cellind,2)),')_',num2str(cell_2),'.tif']);
            multiwaitbar(1,[cellind/size(cellcount,1)],{['Selected (',num2str(round(cellind/size(cellcount,1)*100)),'%): ',...
                num2str(cell_2),' (Mitok), ',num2str(cell_1),' (Perox).']},H)
        end
    end
end
save('cellselected.mat','cellselected');
delete(H.figure)
%% repair cellselected.mat
clc,
mifiles = dir('rawtif_selected_mi\*.tif');
pefiles = dir('rawtif_selected_pe\*.tif');
cellselected = [];
for ii = 1:1:length(mifiles)
    mifilename = mifiles(ii).name;
    splitted = split(string(mifilename),["(",")","_","."]);
    gp = char(splitted(1));
    gp = str2double(gp(end-4:end));
    sb = str2double(splitted(2));
    pg = str2double(splitted(end-1));

    cell_to_select = [gp,sb,pg,0];

    for jj = 1:1:length(pefiles)
        pefilename = pefiles(jj).name;
        splitted2 = split(string(pefilename),["(",")","_","."]);
        gp2 = char(splitted2(1));
        gp2 = str2double(gp2(end-4:end));
        sb2 = str2double(splitted2(2));
        pg2 = str2double(splitted2(end-1));
        if (gp==gp2) && (sb==sb2) && ((pg==pg2+1) || (pg==pg2-1))
            cell_to_select(end) = pg2;
            break
        elseif jj == length(pefiles)
                disp([num2str([ii,jj]),'No matched pefilename found!'])
        end
    end

    cellselected = cat(1,cellselected,cell_to_select);
end
%%
uiprogressdlg