function [imstacksrf,thresholdout] = getnusurface(imstack,rowc,colc,pagec)
%SLICEVIEWER2 Summary of this function goes here
%   Detailed explanation goes here
% imstack = rand(100,100,61);

fig = uifigure("Name","Parallel Stacks Viewer","Position",[800, 300, 640, 400]);
ax1 = uiaxes('Parent',fig,'Units','pixels','Position', [0, 80, 300, 300]);axis(ax1,'tight');
ax2 = uiaxes('Parent',fig,'Units','pixels','Position', [320, 80, 300, 300]);axis(ax2,'tight');
imstack = single(imstack)/max(single(imstack(:)));
ztot = size(imstack,3);
thresholdout = 0.15;
bwstack = imstack>0.15;
imstacksrf = calcsrf(bwstack);
halfz = round(size(bwstack,3)/2+1);
imshow(imstack(:,:,halfz),"Parent",ax1)
imshow(uint8(cat(3,bwstack(:,:,halfz)>0.15,imstacksrf(:,:,halfz),imstacksrf(:,:,halfz))*180),"Parent",ax2)
thres = uieditfield(fig,"numeric","Position",[50,25,100,24],"Limits",[0 1],"Value",0.15,"ValueChangedFcn",@updateim2);
zind = uicontrol(fig,"Style","slider",'SliderStep', [1/ztot, 1/ztot],'Min', 1, 'Max', ztot, 'Value', round(ztot/2+1),"Position",[10,90,620,20],"Callback",@shownow);
uibutton(fig,"Position",[200,25,100,24],"Text","Centroid Connect","ButtonPushedFcn",@getsrf);

waitfor(fig);

    function shownow(hObject,eventdata)
        slider_value = get(hObject,'Value');
        imshow(imstack(:,:,round(slider_value)),"Parent",ax1);title(ax1,round(slider_value) + "/"+ztot);
        imshow(uint8(cat(3,bwstack(:,:,round(slider_value)),...
            imstacksrf(:,:,round(slider_value)),imstacksrf(:,:,round(slider_value)))*250)...
            ,"Parent",ax2);title(ax2,round(slider_value) + "/"+ztot);
    end
    function updateim2(src,event)
        thresholdout = thres.Value;
        disp("Set threshold: "+thres.Value)
        bwstack = imstack>thres.Value;
        imstacksrf = calcsrf(imstack>thres.Value);
        imshow(uint8(cat(3,bwstack(:,:,round(zind.Value)),...
            imstacksrf(:,:,round(zind.Value)),imstacksrf(:,:,round(zind.Value)))*250),"Parent",ax2);
    end
    function getsrf(src,event)
        disp("Set threshold: "+thres.Value)
        thresholdout = thres.Value;
        imstacksrf = centconn(imstacksrf,rowc,colc,pagec);
        imshow(uint8(cat(3,bwstack(:,:,round(zind.Value)),...
            imstacksrf(:,:,round(zind.Value)),imstacksrf(:,:,round(zind.Value)))*250),"Parent",ax2);
        disp("Centroid connection done!")
    end
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