function [Rx,Ry,Rz] = getRxyz(imstack,xythres,mode)
%GETRXYZ Summary of this function goes here
%   Detailed explanation goes here
% thres = 0.5;
imstack = rescale(single(imstack),0,1);
imstack_proj = rescale(max(imstack,[],3),0,1);
[X,Y,Z] = meshgrid(1:1:size(imstack,2),1:1:size(imstack,1),1:1:size(imstack,3));
centroid = [sum(X.*imstack,"all")/sum(imstack,"all"),...
            sum(Y.*imstack,"all")/sum(imstack,"all"),...
            sum(Z.*imstack,"all")/sum(imstack,"all")];
centroid = round(centroid);
% centroid = round(regionprops3(imstack>=0.5,"Centroid").(1)(1,:));
disp(centroid)
imstackxy = imstack(:,:,centroid(3));
imstackxz = squeeze(imstack(centroid(2),:,:))';
% imstackxz = deconvlucy(imstackxz,fspecial("gaussian",size(imstackxz,2),2));
imstackxz = rescale(imstackxz,0,1);
imstackyz = squeeze(imstack(:,centroid(1),:))';
% imstackyz = deconvlucy(imstackyz,fspecial("gaussian",size(imstackxz,2),2));
imstackyz = rescale(imstackyz,0,1);

fig4 = figure(4);clf;
% ax4 = axes(fig4);
subax1 = subplot(2,2,1,Parent=fig4);
imshow(imstackxy,Parent=subax1),title("imstack_proj","Interpreter","none",Parent=subax1);drawnow;
subax2 = subplot(2,2,2,Parent=fig4);
imshow(imstackxz,Parent=subax2),title("imstack_xy","Interpreter","none",Parent=subax2);drawnow;
subax3 = subplot(2,2,3,Parent=fig4);
imshow(imstackyz,Parent=subax3),title("imstack_xz","Interpreter","none",Parent=subax3);drawnow;
subax4 = subplot(2,2,4,Parent=fig4);
imshow(imstack_proj*255,hot(256),Parent=subax4),title("imstack_yz","Interpreter","none",Parent=subax4);drawnow;

% ========== find circles ============
[centers_xy,radii_xy,metric_xy] = imfindcircles(imstackxy,[1,35],Method="PhaseCode");
[centers_xy,radii_xy,metric_xy] = pickupmax(centers_xy,radii_xy,metric_xy);
disp("radii_xy" + radii_xy)
if radii_xy<=10
    % >>>>>>>>>>>>>> Do 3D Gauss fitting <<<<<<<<<<<<<<<<<<<
    if radii_xy>=xythres
        if mode=="wf"
            imstack = deconvlucy(imstack,fspecial3("gaussian",15,[4,4,7]),10);
        end
        [PosRaw,SizeRaw,~] = gaussfit(imstack,'3D',[1 1 1]);
        SizeRaw = SizeRaw/2.35482*3;
        Rx = mean([SizeRaw,15-radii_xy],"all")/2;%SizeRaw(1)/2;
        Ry = mean([SizeRaw,15-radii_xy],"all")/2;%SizeRaw(2)/2;
        Rz = mean([SizeRaw,SizeRaw],"all")/2;%SizeRaw(3)/2;
    else
        if mode=="wf"
            imstack = deconvlucy(imstack,fspecial3("gaussian",15,[4,4,7]),10);
        end
        [PosRaw,SizeRaw,~] = gaussfit(imstack,'3D',[1 1 1]);
        % Rx = mean([SizeRaw,10-radii_xy],"all")/2;%SizeRaw(1)/2;
        % Ry = mean([SizeRaw,10-radii_xy],"all")/2;%SizeRaw(2)/2;
        % Rz = mean(SizeRaw,"all")/2;%SizeRaw(3)/2;
        Rx = SizeRaw(1)/2;
        Ry = SizeRaw(2)/2;
        Rz = SizeRaw(3)/2;
    end
    disp([PosRaw;[Rx,Ry,Rz]]')
    viscircles(subax1,PosRaw(1:2),mean([Rx,Ry]));drawnow;
    viscircles(subax2,PosRaw(1:2:3),mean([Rx,Rz]));drawnow;
    viscircles(subax3,PosRaw(2:3),mean([Ry,Rz]));drawnow;
else
    % >>>>>>>>>>>>>> Do Hough Circle Searching <<<<<<<<<<<<<<<<<<<
    radiirange = 5;
    % [a,b,D] = getaxislength(imstackxy);
    [centers_xz,radii_xz,metric_xz] = imfindcircles(imstackxz,...
                                                    [max(round(radii_xy)-radiirange,1),round(radii_xy)+radiirange],...
                                                    Method="PhaseCode",Sensitivity=0.95);
    [centers_xz,radii_xz,metric_xz] = pickupmax(centers_xz,radii_xz,metric_xz);
    % [a,b,D] = getaxislength(imstackxz);
    [centers_yz,radii_yz,metric_yz] = imfindcircles(imstackyz,...
                                                    [max(round(radii_xy)-radiirange,1),round(radii_xy)+radiirange], ...
                                                    Method="PhaseCode",Sensitivity=0.95);
    [centers_yz,radii_yz,metric_yz] = pickupmax(centers_yz,radii_yz,metric_yz);
    % [a,b,D] = getaxislength(imstackyz);
    disp([centers_xy,radii_xy,metric_xy])
    disp([centers_xz,radii_xz,metric_xz])
    disp([centers_yz,radii_yz,metric_yz])
    viscircles(subax1,centers_xy,radii_xy);drawnow;
    viscircles(subax2,centers_xz,radii_xz);drawnow;
    viscircles(subax3,centers_yz,radii_yz);drawnow;
    % ========== end of find circles ============
    % Rxy = [radii_xy,radii_xy];
    % Rxz = [radii_xz,radii_xz];
    % Ryz = [radii_yz,radii_yz];
    Rx = max([radii_xy,radii_xz]);%(Rxy(1) + Rxz(1))/2;
    Ry = max([radii_xy,radii_yz]);%(Rxy(2) + Ryz(1))/2;
    Rz = max([radii_yz,radii_xz]);%(Rxz(2) + Ryz(2))/2;
end
% ========== find ellipse ============
% [a_xy,b_xy,D_xy] = 
% ========== end of find ellipse ============
Rx = double(Rx);
Ry = double(Ry);
Rz = double(Rz);
end

function [centers,radii,metric] = pickupmax(centers,radii,metric)
if size(centers,1)>1
    centers = centers(radii==max(radii(:)),:);
    radii = radii(radii==max(radii(:)));
    metric = metric(radii==max(radii(:)));
end
end

% function [a,b,D] = getaxislength(im)
% im = rescale(single(im),0,1);
% [~,Y] = meshgrid(1:1:size(im,2),1:1:size(im,1));
% % XC = sum(X.*im,"all")/sum(im,"all");
% YC = round(sum(Y.*im,"all")/sum(im,"all"));
% % centroid = round([XC,YC]);
% impart = zeros(size(im,1),size(im,2),"single");
% windowsize = 8;
% YC = min(max(windowsize+1,YC),size(im,1)-windowsize);
% impart(YC-windowsize:YC+windowsize,:) = im(YC-windowsize:YC+windowsize,:);
% impart = rescale(impart,0,1);
% impartedge = edge(impart.*single(impart>0.5),"canny");
% for ii = YC-windowsize:YC+windowsize
%     k = find(impartedge(ii,:)>0);
%     if length(k)>2
%         for kk = 2:length(k)-1
%             impartedge(ii,k(kk))=0;
%         end
%     end
% end
% % fig = figure;ax = axes(fig);
% % imshowpair(impart,impartedge,"montage",Parent=ax);
% a = 10;
% b = 10;
% D = 10;
% end