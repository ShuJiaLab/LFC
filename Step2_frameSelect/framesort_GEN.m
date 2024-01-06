clear,clc
load("cellcount.mat","cellcount")
file_prefix = 'JK_minu60_';
%% make folders
mkdir('./z_unused_000')
mkdir('./z_single_nucleus')
mkdir('./z_multi_nucleus')
mkdir('./z_multi_mito')
%% Sort single frames
% These frames may contain nucleus only.
% ------------------------------------------
for ii = 1:1:length(cellcount)
    filename = ['./rawtif_selected/',file_prefix,num2str(cellcount(ii,1),'%05d'),...
            '(',num2str(cellcount(ii,2)),')','_',num2str(cellcount(ii,3)),'.tif'];
    if cellcount(ii,5)==1
        movefile(filename,'./z_single_nucleus');
    end
end
disp([length(dir('./rawtif_selected/'))-2,length(dir('./z_single_nucleus/'))-2])
% Now 'z_single_nucleus' folder has single nucleus, or single mitochondria
% or just cell residues/dusts, which require further selection.
%% Sort multiple frmes
clc,
try
    load("loginfo.mat","loginfo");
catch
    loginfo.cellsortedpt = 0;
    loginfo.cellreconpt = 0;
    save("loginfo.mat","loginfo");
end
% ------------------------------------------
H = waitbar(0,num2str([length(dir('./rawtif_selected/'))-2,...
    length(dir('./z_multi_nucleus/'))-2,length(dir('./z_multi_mito/'))-2]));
for ii = loginfo.cellsortedpt+1:1:length(cellcount)   
    if cellcount(ii,5)~=1
        imstack = [];
        waitbar(ii/length(cellcount),H,['Getting image batch: ',num2str(cellcount(ii,:))]);
        for jj = cellcount(ii,3):1:cellcount(ii,4)
            filename = ['./rawtif_selected/',file_prefix,num2str(cellcount(ii,1),'%05d'),...
            '(',num2str(cellcount(ii,2)),')','_',num2str(jj),'.tif'];
            imstack = cat(3,imstack,imread(filename));
        end
        colsize = 4;
        imbatch = stack2batch(imstack,colsize);
        [row,col,NumOfFrames] = size(imstack);
        while 1
            try
                imorder = batch2order(imbatch,row,col,colsize);
            catch
                delete(H)
                return
            end
            if length(imorder)==1
                filename = ['./rawtif_selected/',file_prefix,num2str(cellcount(ii,1),'%05d'),...
                    '(',num2str(cellcount(ii,2)),')','_',num2str(cellcount(ii,3)+imorder-1),'.tif'];
                movefile(filename,'./z_single_nucleus');
                disp({cellcount(ii,3)+imorder-1,'./single_nucleus/'})
            elseif length(imorder)==2
                filename = ['./rawtif_selected/',file_prefix,num2str(cellcount(ii,1),'%05d'),...
                    '(',num2str(cellcount(ii,2)),')','_',num2str(cellcount(ii,3)+imorder(1)-1),'.tif'];
                movefile(filename,'./z_multi_nucleus');
                filename = ['./rawtif_selected/',file_prefix,num2str(cellcount(ii,1),'%05d'),...
                    '(',num2str(cellcount(ii,2)),')','_',num2str(cellcount(ii,3)+imorder(2)-1),'.tif'];
                movefile(filename,'./z_multi_mito');
                disp({cellcount(ii,3)+imorder(1)-1,'./multi_nucleus/';...
                      cellcount(ii,3)+imorder(2)-1,'./multi_mito/'})
            else
                disp({cellcount(ii,3),cellcount(ii,4),'[No cell selected!]'})
            end
            answer = questdlg('Multiple cells in the frame series?','Multi-choices','Yes','No','No');
            if strcmp(answer,'No')
                break
            end
        end
    end

    waitbar(ii/length(cellcount),H,num2str([length(dir('./rawtif_selected/'))-2,...
        length(dir('./z_multi_nucleus/'))-2,length(dir('./z_multi_mito/'))-2]))

    loginfo.cellsortedpt = ii;
    loginfo.cellreconpt = 0;
    save("loginfo.mat","loginfo");
end
delete(H)


%% ======================= Utilities =======================
% Here contains all the subfunctions that are used in the main function.
% ==========================================================
function imbatch = stack2batch(imstack,colsize)
[row,col,NumOfFrames] = size(imstack);
rowsize = ceil(NumOfFrames/colsize);
imbatch = zeros(row*rowsize,col*colsize);
imstack = single(imstack)./max(single(imstack),[],[1,2]);
for ii = 1:1:NumOfFrames
    row2put = ceil(ii/colsize);
    col2put = round(ii - (row2put-1)*colsize);
    imbatch((row2put-1)*row+1:row2put*row,...
        (col2put-1)*col+1:col2put*col) = imstack(:,:,ii);
end
end
% --------------------------------
function imorder = batch2order(imbatch,row,col,colsize)
figure(99),imshow(imbatch);
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