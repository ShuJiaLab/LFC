function [nustackaligned,mistackaligned,...
          rowc,colc,pagec,tomove,...
          nustacksrf,mistacksrf,mistack_nu,mistack_all] = NuMiCore(nustack,mistack)
%NUMICORE Summary of this function goes here
%   Detailed explanation goes here
% nustack = rand(100,100,10);
% mistack = rand(100,100,10);
% -------------------------------------------------
nustack = single(nustack)/max(single(nustack(:)));
mistack = single(mistack)/max(single(mistack(:)));
nustackaligned = nustack;
mistackaligned = mistack;
nustackaligned2 = nustackaligned;
mistackaligned2 = mistackaligned;
nustacksrf = nustackaligned;
mistacksrf = mistackaligned;
mistack_nu = nustackaligned;
mistack_all = mistackaligned;
fig = uifigure("Name","Parallel Stacks Viewer","Position",[300, 200, 1280, 700]);
% ==============================================================================================================
ax1 = uiaxes('Parent',fig,'Units','pixels','Position', [0, 390, 300, 300]);axis(ax1,'tight');
ax2 = uiaxes('Parent',fig,'Units','pixels','Position', [320, 390, 300, 300]);axis(ax2,'tight');
ztot = size(nustack,3);
halfz = round(ztot/2+1);
rowc=1;colc=1;pagec=halfz;
mistack = single(mistack)/max(single(mistack(:)));
imshow(nustack(:,:,halfz),"Parent",ax1);title(ax1,round(halfz) + "/"+ztot);
imshow(mistack(:,:,halfz),"Parent",ax2);title(ax2,round(halfz) + "/"+ztot);
zind = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot, 1/ztot],...
    'Min', 1, 'Max', ztot, 'Value', round(ztot/2+1),"Position",[10,380,620,20],...
    "Callback",@(zind,eventdata)shownow(zind,nustack,mistack));
% ==============================================================================================================
zprofile_max_1 = rescale(squeeze(max(nustack,[],[1,2])),0,1);
zprofile_sum_1 = rescale(squeeze(sum(nustack,[1,2])),0,1);
zprofile_max_2 = rescale(squeeze(max(mistack,[],[1,2])),0,1);
zprofile_sum_2 = rescale(squeeze(sum(mistack,[1,2])),0,1);
tomove = 0;
z0 = halfz;
ax3 = uiaxes('Parent',fig,'Units','pixels','Position', [640, 410, 300, 270]);axis(ax3,'tight');
plot(1:1:length(zprofile_max_1),zprofile_max_1,"r-",...
     1:1:length(zprofile_sum_1),zprofile_sum_1,"b-",...
     1:1:length(zprofile_max_2),zprofile_max_2,"r--",...
     1:1:length(zprofile_sum_2),zprofile_sum_2,"b--","Parent",ax3);
focalplane = uieditfield(fig,"numeric","Position",[670,380,100,24],"Value",z0,"ValueChangedFcn",@selectz);
tomovebox = uicheckbox(fig,"Text","Move to trash?","Position",[820,380,100,24],"Value",0,"ValueChangedFcn",@movechanged);
% ==============================================================================================================
nustackproj = mean(single(nustack),3);
mistackproj = mean(single(mistack),3);
ax4 = uiaxes('Parent',fig,'Units','pixels','Position', [960, 390, 300, 300]);axis(ax4,'tight');
imshowpair(nustackproj,mistackproj,"Parent",ax4);
goupdown = uieditfield(fig,"numeric","Position",[980,380,100,24],"Value",0,"ValueChangedFcn",@moveimage);
goleftright = uieditfield(fig,"numeric","Position",[1140,380,100,24],"Value",0,"ValueChangedFcn",@moveimage);
% ==============================================================================================================
ax5 = uiaxes('Parent',fig,'Units','pixels','Position', [0, 70, 300, 300]);axis(ax5,'tight');
ax6 = uiaxes('Parent',fig,'Units','pixels','Position', [320, 70, 300, 300]);axis(ax6,'tight');
thresh= 0.15; bwstacknu = nustackaligned2>thresh;
nustacksrf = calcsrf(bwstacknu);
ztot2 = size(bwstacknu,3);
halfz2 = round(ztot2/2+1);
imshow(nustackaligned2(:,:,halfz2),"Parent",ax5);title(ax5,halfz2 + "/"+ztot2);
imshow(uint8(cat(3,bwstacknu(:,:,halfz2)>thresh,nustacksrf(:,:,halfz2),bwstacknu(:,:,halfz2)>thresh)*180),"Parent",ax6);
title(ax6,halfz2 + "/"+ztot2);
thres = uieditfield(fig,"numeric","Position",[50,20,100,24],"Limits",[0 1],"Value",thresh,"ValueChangedFcn",@updateim2);
zind2 = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot2, 1/ztot2],...
    'Min', 1, 'Max', ztot2, 'Value', round(ztot2/2+1),"Position",[10,60,620,20],"Callback",@shownow2);
uibutton(fig,"Position",[200,20,120,24],"Text","Centroid Connect","ButtonPushedFcn",@getsrf);
% ==============================================================================================================
ax7 = uiaxes('Parent',fig,'Units','pixels','Position', [640, 70, 300, 300]);axis(ax7,'tight');
ax8 = uiaxes('Parent',fig,'Units','pixels','Position', [960, 70, 300, 300]);axis(ax8,'tight');
thresh2= 0.15; bwstackmi = mistackaligned2>thresh2;
mistacksrf = bwstackmi;
ztot3 = size(bwstackmi,3);
halfz3 = round(ztot3/2+1);
imshow(mistackaligned2(:,:,halfz3),"Parent",ax7);title(ax7,halfz3 + "/"+ztot3);
imshow(uint8(cat(3,bwstackmi(:,:,halfz3)>thresh2,mistacksrf(:,:,halfz3),bwstackmi(:,:,halfz3)>thresh2)*180),"Parent",ax8);
title(ax8,halfz3 + "/"+ztot3);
thres2 = uieditfield(fig,"numeric","Position",[700,20,100,24],"Limits",[0 1],"Value",thresh2,"ValueChangedFcn",@updateim3);
zind3 = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot3, 1/ztot3],...
    'Min', 1, 'Max', ztot3, 'Value', round(ztot3/2+1),"Position",[650,60,620,20],"Callback",@shownow3);

waitfor(fig)
nustackaligned = nustackaligned2;
mistackaligned = mistackaligned2;
% disp([size(mistackaligned2);size(nustacksrf)])
mistack_nu = mistackaligned2.*double(nustacksrf)*60000;
mistack_all = mistackaligned2.*double(mistacksrf)*60000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function shownow(hObject,imstack1,imstack2)
        slider_value = get(hObject,'Value');
        imshow(imstack1(:,:,round(slider_value)),"Parent",ax1);title(ax1,round(slider_value) + "/"+ztot);
        imshow(imstack2(:,:,round(slider_value)),"Parent",ax2);title(ax2,round(slider_value) + "/"+ztot);
    end

    function movechanged(src,event)
        tomove = tomovebox.Value;
        if tomove
            rowc=0;colc=0;pagec=0;
            focalplane.Value = 0;
            disp("<Move to Trash> Set centroid: ["+rowc+", "+colc+", "+pagec+"]")
        else
            pagec=halfz;
            focalplane.Value = halfz;
            disp("<Restore from Trash> Set centroid: ["+rowc+", "+colc+", "+pagec+"]")
        end
    end

    function selectz(src,event)
        z0 = focalplane.Value;pagec=z0;
        disp("Centroid: ["+rowc+", "+colc+", "+pagec+"]")
    end

    function moveimage(src,event)
        dy = goupdown.Value;
        dx = goleftright.Value;
        nustackaligned = circshift(nustack,[dy,dx]);
        nustackproj = mean(single(nustackaligned),3);
        imshowpair(nustackproj,mistackproj,"Parent",ax4);drawnow;
        [rowc,colc,~] = getcent(nustackaligned,0.1);pagec = z0;
        disp("Centroid: ["+rowc+", "+colc+", "+pagec+"]")
        nustackaligned2 = nustackaligned(:,:,pagec-15:pagec+15);
        mistackaligned2 = mistackaligned(:,:,pagec-15:pagec+15);
        ztot2 = size(nustackaligned2,3);
        halfz2 = round(ztot2/2+1);
        bwstacknu = nustackaligned2>thresh;
        nustacksrf = calcsrf(bwstacknu);
        zind2 = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot2, 1/ztot2],...
            'Min', 1, 'Max', ztot2, 'Value', round(ztot2/2+1),...
            "Position",[10,60,620,20],"Callback",@shownow2);
        ztot3 = size(mistackaligned2,3);
        halfz3 = round(ztot3/2+1);
        bwstackmi = mistackaligned2>thresh2;
        mistacksrf = bwstackmi;
        zind3 = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot3, 1/ztot3],...
            'Min', 1, 'Max', ztot3, 'Value', round(ztot3/2+1),...
            "Position",[650,60,620,20],"Callback",@shownow3);
    end

    function shownow2(hObject,event)
        slider_value = get(hObject,'Value');
        imshow(nustackaligned2(:,:,round(slider_value)),"Parent",ax5);title(ax5,round(slider_value) + "/"+ztot2);
        imshow(uint8(cat(3,bwstacknu(:,:,round(slider_value)),...
            nustacksrf(:,:,round(slider_value)),bwstacknu(:,:,round(slider_value)))*250)...
            ,"Parent",ax6);title(ax6,round(slider_value) + "/"+ztot2);
    end

    function updateim2(src,event)
        thresh = get(src,'Value');
        disp("Set threshold: "+thresh)
        bwstacknu = nustackaligned2>thresh;
        nustacksrf = calcsrf(bwstacknu);
        imshow(uint8(cat(3,bwstacknu(:,:,round(zind2.Value)),...
            nustacksrf(:,:,round(zind2.Value)),bwstacknu(:,:,round(zind2.Value)))*250),"Parent",ax6);
    end

    function getsrf(src,event)
        disp("Set threshold: "+thres.Value)
        thresh = thres.Value;
        nustacksrf = centconn(nustacksrf,rowc,colc,halfz2);
        imshow(uint8(cat(3,bwstacknu(:,:,round(zind2.Value)),...
            nustacksrf(:,:,round(zind2.Value)),bwstacknu(:,:,round(zind2.Value)))*250),"Parent",ax6);
        disp("Centroid connection done!")
    end

     function shownow3(hObject,event)
        slider_value = get(hObject,'Value');
        imshow(mistackaligned2(:,:,round(slider_value)),"Parent",ax7);title(ax7,round(slider_value) + "/"+ztot3);
        imshow(uint8(cat(3,bwstackmi(:,:,round(slider_value)),...
            mistacksrf(:,:,round(slider_value)),bwstackmi(:,:,round(slider_value)))*250)...
            ,"Parent",ax8);title(ax8,round(slider_value) + "/"+ztot3);
    end

    function updateim3(src,event)
        thresh2 = get(src,'Value');
        disp("Set threshold: "+thresh2)
        bwstackmi = mistackaligned2>thresh2;
        mistacksrf = bwstackmi;
        imshow(uint8(cat(3,bwstackmi(:,:,round(zind3.Value)),...
            mistacksrf(:,:,round(zind3.Value)),bwstackmi(:,:,round(zind3.Value)))*250),"Parent",ax8);
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function bwstacksrf = calcsrf(bw)
strelsize = 3;
se135 = strel('line',strelsize,135);
se90 = strel('line',strelsize,90);
se45 = strel('line',strelsize,45);
se0 = strel('line',strelsize,0);
BWsdil = imdilate(bw,[se135 se90 se45 se0]);
bwstacksrf = imclose(BWsdil,strel('disk',20));
end

function imstacksrfout = centconn(imstacksrf,rowc,colc,pagec)
[rows,cols,pages] = ind2sub(size(imstacksrf),find(imstacksrf>0));
imstacksrfout = imstacksrf;
for ii = 1:1:length(rows)
    if rows(ii)~=rowc
        m_n = (cols(ii)-colc)/(rows(ii)-rowc);
        p_n = (pages(ii)-pagec)/(rows(ii)-rowc);
        xi = m_n*(1:1:((rows(ii)-rowc)))+colc;
        yi = (1:1:(rows(ii)-rowc)) + rowc;
        zi = p_n*(1:1:((rows(ii)-rowc)))+pagec;
        xi = round(xi);yi = round(yi);zi = round(zi);
        for jj = 1:1:length(xi)
            imstacksrfout(yi(jj),xi(jj),zi(jj)) = 1;
        end
    else
        if cols(ii)~=colc
            p_m = (pages(ii)-pagec)/(cols(ii)-colc);
            zi = p_m*(1:1:((cols(ii)-colc)))+pagec;
            xi = (1:1:(cols(ii)-colc))+colc;
            xi = round(xi);yi = round(yi);zi = round(zi);
            for jj = 1:1:length(xi)
                imstacksrfout(rowc,xi(jj),zi(jj)) = 1;
            end
        else
            zi = 1:1:((pages(ii)-pagec))+pagec;
            xi = round(xi);yi = round(yi);zi = round(zi);
            xi = round(xi);yi = round(yi);zi = round(zi);
            for jj = 1:1:length(xi)
                imstacksrfout(rowc,colc,zi(jj)) = 1;
            end
        end

    end
end
end