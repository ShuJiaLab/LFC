clear,clc

im1name = ['rawtif_selected_r\bead_flf_r00001(1)_','827','.tif'];
im2name = ['rawtif_selected_r\bead_flf_r00001(1)_','853','.tif'];

im1 = imread(im1name);
im2 = imread(im2name);

bkg = 90;
thres_2 = 150;
g1 = [-1,0,1;-2,0,2;-1,0,1]; % Sobel Gx kernel
g2 = -g1';

im1_cp = (double(im1)-bkg)/max(double(im1(:)));
im1_cp(im1_cp<=0) = 0;
im2_cp = (double(im2)-bkg)/max(double(im2(:)));
im2_cp(im2_cp<=0) = 0;
t1 = thres_2/double(max(im1(:)));
t2 = thres_2/double(max(im2(:)));
t = min([max([t1,t2]),0.5]);
% t = 0.65;

se = strel('disk',10);
im1_bw = imbinarize(im1_cp,t);
im1_closed = imclose(im1_bw,se);
sum1 = sum(abs(im1_closed-im1_bw),'all')/sum(im1_closed,'all');
% f1 = sum((filter2(im1_cp,[0 1 0;1 -4 1;0 1 0])).^2,'all');

im2_bw = imbinarize(im2_cp,t);
im2_closed = imclose(im2_bw,se);
sum2 = sum(abs(im2_closed-im2_bw),'all')/sum(im2_closed,'all');
% f2 = sum((filter2(im2_cp,[0 1 0;1 -4 1;0 1 0])).^2,'all');

figure(1),
subplot(321),imshow(im1_cp),title(num2str(sum1));
subplot(322),imshow(im2_cp),title(num2str(sum2));
subplot(323),imshow(im1_bw);
subplot(324),imshow(im2_bw);
subplot(325),imshow(im1_closed);
subplot(326),imshow(im2_closed);
%%
scrfiles = dir('.\rawtif_selected_b\*.tif');
bkg = 90;
thres_2 = 150;
h = waitbar(0,'Processing...');
for ii = 1:length(scrfiles)
    im1 = imread(['.\rawtif_selected_b\',scrfiles(ii).name]);
    im1_cp = (double(im1)-bkg)/max(double(im1(:)));
    im1_cp(im1_cp<=0) = 0;
    t1 = thres_2/double(max(im1(:)));
    t = min([t1,0.5]);
    se = strel('disk',10);
    im1_bw = imbinarize(im1_cp,t);
    im1_closed = imclose(im1_bw,se);
    sum1 = sum(abs(im1_closed-im1_bw),'all')/sum(im1_closed,'all');
    if sum1>0.15
%         movefile(['.\rawtif_selected_b_focused\',scrfiles(ii).name],'.\rawtif_selected_b_defocused\')
        imwrite(im1,['.\rawtif_selected_b_defocused\',scrfiles(ii).name])
        removed = ' (x) ';
    else
        removed = ' ( ) ';
        imwrite(im1,['.\rawtif_selected_b_focused\',scrfiles(ii).name])
    end
    waitbar(ii/length(scrfiles),h,...
        [scrfiles(ii).name,removed,...
        num2str(ii),'/',num2str(length(scrfiles)),', (',num2str(round(ii/length(scrfiles)*100,2)),'%)'])
end
delete(h)
disp('Done!')
%% =============================================
% ==============================================
scrfiles = dir('.\rawtif_selected_r_focused\*.tif');
frameres = [];
h = waitbar(0,'Processing...');
for ii = 1:length(scrfiles)
    scrfilename = split(scrfiles(ii).name,'_');
    scrfilename3 = split(scrfilename{3},["(",")"]);
    scrfilename4 = split(scrfilename{4},'.');
    fileind = str2double(scrfilename3{1}(2:end));
    tifnum = str2double(scrfilename3{2});
    pageind = str2double(scrfilename4{1});
    frameres = cat(1,frameres,[fileind,tifnum,pageind]);
    waitbar(ii/length(scrfiles),h,[num2str(round(ii/length(scrfiles)*100,2)),'%']);
end
delete(h)
frameres = sortrows(frameres);
%%
clc,
cellcount = [];
H = waitbar(0,'cellcount: []');
for ind = 1:1:(size(frameres,1)-1)
    if ind == 1
        ctn = (frameres(ind+1,3) == frameres(ind,3)+1);
        cellcount = [frameres(ind,1),frameres(ind,2),frameres(ind,3),frameres(ind,3),1];
        if ~ctn
            cellcount = cat(1,cellcount,...
                [frameres(ind+1,1),frameres(ind+1,2),frameres(ind+1,3),frameres(ind+1,3),1]);
        end
        waitbar(ind/(size(frameres,1)-1),H,['#',num2str(ind),': ',num2str(cellcount(end,:))]);
    else
        ctn = (frameres(ind+1,3) == (frameres(ind,3)+1));
        if ~ctn
            cellcount(end,4) = frameres(ind,3);
            cellcount(end,5) = cellcount(end,4) - cellcount(end,3) + 1;
            cellcount = cat(1,cellcount,...
                [frameres(ind+1,1),frameres(ind+1,2),frameres(ind+1,3),frameres(ind+1,3),1]);
            waitbar(ind/(size(frameres,1)-1),H,['#',num2str(ind),': ',num2str(cellcount(end-1,:))]);
        else
            if ind == (size(frameres,1)-1)
            cellcount(end,5) = cellcount(end,4) - cellcount(end,3) + 1;
            end
        end
    end
end
% cellcount(end,4) = frameRes2(ind,3);
delete(H)
% cellcount(end,5) = cellcount(end,4) - cellcount(end,3) + 1;
%%
clc,
frameRes3 = find(cellcount(:,5)>1);
h = waitbar(0,'Processing...');
toremove = [];
for ii = 1:length(cellcount)
    cents = [];
    for jj = cellcount(ii,3):1:cellcount(ii,4)
        im = imread(['.\rawtif_selected_r_focused\bead_flf_r',num2str(cellcount(ii,1),'%05d'),'(',num2str(cellcount(ii,2)),')_',num2str(jj),'.tif']);
        im = im(501:900,211:810);
        im = (double(im)-90)/(max(double(im(:)))-90);
        im(im<=0) = 0;
        im = imbinarize(im,0.5);
        [X,Y] = meshgrid(1:1:size(im,2),1:1:size(im,1));
        centv = sqrt((sum(im.*X,'all')/sum(im,'all')-size(im,2)/2)^2 + (sum(im.*Y,'all')/sum(im,'all')-size(im,1)/2)^2);
        cent = [jj,centv,double(centv>=200)];
        cents = cat(1,cents,cent);
    end
%     disp(cents)
%     disp('=========================')
    [~,centmini] = min(cents(:,2));
    toremove = cat(1,toremove,[cellcount(ii,3),cellcount(ii,4),cents(centmini,1)]);
    if cellcount(ii,3)~=cellcount(ii,4)
        for jj = cellcount(ii,3):1:cellcount(ii,4)
            im_filename = ['.\rawtif_selected_r_focused\bead_flf_r',num2str(cellcount(ii,1),'%05d'),'(',num2str(cellcount(ii,2)),')_',num2str(jj),'.tif'];
            if (jj ~= cents(centmini,1)) || (cents(cents(:,1)==jj,end)==1)
                movefile(im_filename,'.\rawtif_selected_r_decentered\')
                disp(['moved: ',im_filename,' | ',num2str(cents(cents(:,1)==jj,end))])
            end
        end
    else
        im_filename = ['.\rawtif_selected_r_focused\bead_flf_r',num2str(cellcount(ii,1),'%05d'),'(',num2str(cellcount(ii,2)),')_',num2str(cellcount(ii,3)),'.tif'];
        if (cents(cents(:,1)==cellcount(ii,3),end)==1)%(jj ~= cents(centmini,1))||
            movefile(im_filename,'.\rawtif_selected_r_decentered\')
            disp(['moved: ',im_filename,' | ',num2str(cents(cents(:,1)==jj,end))])
        end
    end
%     waitbar(find(frameRes3==ii)/length(frameRes3),h,['Processing... [',num2str(cents(:,1)'),'], ',num2str(cents(centmini,1))])
    waitbar(ii/length(cellcount),h,['Processing... [',num2str(cents(:,1)'),'], ',num2str(cents(centmini,1))])
end
delete(h)
disp('Done!')
%% ===============================================================================================
clc,
working_folder = '.\rawtif_selected_r_focused\';
scrfiles = dir([working_folder,'*.tif']);
working_files = [];
h = waitbar(0,'Getting working file indices...');
for ii = 1:length(scrfiles)
    scrfilename = split(scrfiles(ii).name,'_');
    scrfilename3 = split(scrfilename{3},["(",")"]);
    scrfilename4 = split(scrfilename{4},'.');
    fileind = str2double(scrfilename3{1}(2:end));
    tifnum = str2double(scrfilename3{2});
    pageind = str2double(scrfilename4{1});
    working_files = cat(1,working_files,[fileind,tifnum,pageind]);
    waitbar(ii/length(scrfiles),h,['Getting working file indices... (',num2str(round(ii/length(scrfiles)*100,2)),'%)']);
end
delete(h)
working_files = sortrows(working_files);
%% ************************************************************************************************
clc,
GaussMod = @(BETA,x)...      % P is related to the size of the peak, greater->narrower
            BETA(1)*(...
            exp( - 1./BETA(2).*((x(:,1)-BETA(3)).^2 )) .* ...
            exp( - 1./BETA(4).*((x(:,2)-BETA(5)).^2 )) ...
            )+BETA(6);
mserrs1 = [];
h = waitbar(0,'Performing 2D Gauss Fitting...');
for ii = 1:1:length(working_files)
    filename = [working_folder,'bead_flf_r',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)),'.tif'];
    im = imread(filename);
    
    im = double(im(501:1000,211:810));
    im = (im-90)/(max(im(:))-90);
    im(im<=0) = 0;
%     im = sqrt(im);
    [rows,cols] = find(imbinarize(im,0.75) == 1);
    rowc = round(mean(rows,"all"));
    colc = round(mean(cols,"all"));
    rc = sqrt(sum(imbinarize(im,0.5),"all")/pi)*2;
    try
        im = im((rowc-49):(rowc+50),(colc-49):(colc+50));
        [CoX,CoY] = meshgrid(1:1:size(im,2),1:1:size(im,2));
        X = [CoX(:),CoY(:)];
        Y = im(:);
        StartP = [1,(rc^2)/log(16)*2.35,50,(rc^2)/log(16)*2.35,50,0];
    catch
        im = im(max([1,(rowc-99)]):min([(rowc+100),size(im,1)]),max([1,(colc-99)]):min([(colc+100),size(im,2)]));
        [CoX,CoY] = meshgrid(1:1:size(im,2),1:1:size(im,1));
        X = [CoX(:),CoY(:)];
        Y = im(:);
        StartP = [1,(rc^2)/log(16)*2.35,colc,(rc^2)/log(16)*2.35,rowc,0];
    end
    try
    [beta,~,~,~,MSErr,~] = nlinfit(X,Y,GaussMod,StartP');
    if MSErr>0.001
        movefile(filename,'.\rawtif_selected_r_defocused\')
        disp(['moved: ',filename,' | ',num2str(cents(cents(:,1)==jj,end))])
    end
    catch
        MSErr = -1e-4;
        movefile(filename,'.\rawtif_selected_r_defocused\')
        disp(['moved: ',filename,' | ',num2str(cents(cents(:,1)==jj,end))])
    end
    mserrs1 = cat(1,mserrs1,MSErr);
    disp([filename,'    ',num2str(MSErr)]);
    waitbar(ii/length(working_files),h,['Performing 2D Gauss Fitting... (',num2str(round(ii/length(working_files)*100,2)),'%) MSE: ',num2str(round(MSErr,5))]);
end
delete(h)
disp('Fitting Done!')
%% ==============================================================================
% ===============================================================================
% ===============================================================================
clc,
working_folder = '.\rawtif_selected_b_focused\';
% working_folder_0 = '.\rawtif_selected_b_focused\';
scrfiles = dir([working_folder,'*.tif']);
working_files = [];
h = waitbar(0,'Getting working file indices...');
for ii = 1:length(scrfiles)
    scrfilename = split(scrfiles(ii).name,'_');
    scrfilename3 = split(scrfilename{3},["(",")"]);
    scrfilename4 = split(scrfilename{4},'.');
    fileind = str2double(scrfilename3{1}(2:end));
    tifnum = str2double(scrfilename3{2});
    pageind = str2double(scrfilename4{1});
    working_files = cat(1,working_files,[fileind,tifnum,pageind]);
    waitbar(ii/length(scrfiles),h,['Getting working file indices... (',num2str(round(ii/length(scrfiles)*100,2)),'%)']);
end
delete(h)
working_files = sortrows(working_files);
%% ************************************************************************************************
clc,
GaussMod = @(BETA,x)...      % P is related to the size of the peak, greater->narrower
            BETA(1)*(...
            exp( - 1./BETA(2).*((x(:,1)-BETA(3)).^2 )) .* ...
            exp( - 1./BETA(4).*((x(:,2)-BETA(5)).^2 )) ...
            )+BETA(6);
mserrs1 = [];
h = waitbar(0,'Performing 2D Gauss Fitting...');
for ii = 1:1:length(working_files)
    filename = [working_folder,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)),'.tif'];
    filename_m3 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)-3),'.tif'];
    filename_m2 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)-2),'.tif'];
    filename_m1 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)-1),'.tif'];
    filename_p1 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)+1),'.tif'];
    filename_p2 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)+2),'.tif'];
    filename_p3 = [working_folder_0,'bead_flf_b',num2str(working_files(ii,1),'%05d'),'(',num2str(working_files(ii,2)),')_',num2str(working_files(ii,3)+3),'.tif'];
    im = imread(filename);
    
    im = double(im(501:1000,211:810));
    im = (im-90)/(max(im(:))-90);
    im(im<=0) = 0;
%     im = sqrt(im);
    [rows,cols] = find(imbinarize(im,0.75) == 1);
    rowc = round(mean(rows,"all"));
    colc = round(mean(cols,"all"));
    rc = sqrt(sum(imbinarize(im,0.5),"all")/pi)*2;
    try
        im = im((rowc-49):(rowc+50),(colc-49):(colc+50));
        [CoX,CoY] = meshgrid(1:1:size(im,2),1:1:size(im,2));
        X = [CoX(:),CoY(:)];
        Y = im(:);
        StartP = [1,(rc^2)/log(16)*2.35,50,(rc^2)/log(16)*2.35,50,0];
    catch
        im = im(max([1,(rowc-99)]):min([(rowc+100),size(im,1)]),max([1,(colc-99)]):min([(colc+100),size(im,2)]));
        [CoX,CoY] = meshgrid(1:1:size(im,2),1:1:size(im,1));
        X = [CoX(:),CoY(:)];
        Y = im(:);
        StartP = [1,(rc^2)/log(16)*2.35,colc,(rc^2)/log(16)*2.35,rowc,0];
    end
    try
    [beta,~,~,~,MSErr,~] = nlinfit(X,Y,GaussMod,StartP');
    if MSErr<0.001
        if ~exist(filename_m3,'file') && ~exist(filename_m2,'file') &&...
                ~exist(filename_m1,'file') && ~exist(filename_p1,'file') &&...
                ~exist(filename_p2,'file') && ~exist(filename_p3,'file')
            movefile(filename,'.\rawtif_selected_b_focused\')
            disp(['moved: ',filename])
        end
    end
    catch
        MSErr = -1e-4;
    end
    mserrs1 = cat(1,mserrs1,MSErr);
    disp([filename,'    ',num2str(MSErr)]);
    waitbar(ii/length(working_files),h,['Performing 2D Gauss Fitting... (',num2str(round(ii/length(working_files)*100,2)),'%) MSE: ',num2str(round(MSErr,5))]);
end
delete(h)
disp('Fitting Done!')
%% ^^^^^^^^^^^^^^^^^^^^^^^^ BEAD SORTING ^^^^^^^^^^^^^^^^^^^^^^^^^^^
%  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
clc,
h = waitbar(0,'Sorting beads...');
figure(3),
for ii = 1:1:length(working_files)
    disp('-----------------')
    disp(scrfiles(ii).name)
    imraw = imread([scrfiles(ii).folder,'\',scrfiles(ii).name]);
    try
        imraw = imraw(479:978,262:761);
%         imraw = imresize(imraw,275/125);
        [rows,cols] = find((imraw-90)>0.5*max(imraw(:)-90));
        rowc1 = mean(rows,'all');
        colc1 = mean(cols,'all');
        [~,maxi] = max(imraw(:));
        [rowc,colc] = ind2sub(size(imraw),maxi);
        if abs(colc1-colc)>20 || abs(rowc1-rowc)>20
            colc1 = colc;
            rowc1 = rowc;
        end
        rows2proc1 = round(max(1,rowc1-63):min(rowc1+64,size(imraw,1)));
        cols2proc1 = round(max(1,colc1-63):min(colc1+64,size(imraw,2)));
        imraw = imraw(rows2proc1,cols2proc1);
%         imshow(double(imraw)/max(double(imraw(:)))),
         title([scrfiles(ii).name,' size: -- nm']),
        drawnow;
        [~,beadsizes,~] = gaussfit(imraw,'2D',[145 145]);
        beadsize = mean(beadsizes);
        if beadsize<900
            copyfile([scrfiles(ii).folder,'\',scrfiles(ii).name],...
                'rawtif_selected_b_selected\rawtif_selected_b_200\');
        elseif beadsize<1500
            copyfile([scrfiles(ii).folder,'\',scrfiles(ii).name],...
                'rawtif_selected_b_selected\rawtif_selected_b_1000\');
        elseif beadsize<3000
            copyfile([scrfiles(ii).folder,'\',scrfiles(ii).name],...
                'rawtif_selected_b_selected\rawtif_selected_b_2000\');
        else
            copyfile([scrfiles(ii).folder,'\',scrfiles(ii).name],...
                'rawtif_selected_b_selected\rawtif_selected_b_4000\');
        end
        scatter(beadsize,log(double(max(imraw(:)))),12,[0.5 0.5 1],'filled')
        title([scrfiles(ii).name,' size: ',num2str(round(beadsize)),' nm']),
        drawnow,
        hold on
    catch ME
        copyfile([scrfiles(ii).folder,'\',scrfiles(ii).name],...
                'rawtif_selected_b_selected\rawtif_selected_b_nan\');
        scatter(beadsize,log(double(max(imraw(:)))),12,[1 0.5 0.5],'o')
        title([scrfiles(ii).name,' not sorted!']),
        drawnow,
        hold on
        disp([scrfiles(ii).name,' - Not sorted!'])
%         rethrow(ME)
    end
    waitbar(ii/length(working_files),h,['Sorting beads... (',...
        num2str(round(ii/length(working_files),2)),')']);
end
hold off
delete(h)