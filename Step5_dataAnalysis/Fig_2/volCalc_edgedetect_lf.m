clear,clc
warning off
% flfm_folder = '.\rawtif_selected_r_selected\rawtif_selected_r_200_selected_offsetsub_acsn_ccut\';
recon_folder = '.\PSFFLFint_20230706_Red_Gly_10um-bead200lf_ccut_single\';
% recon_folder = '.\bead200_single\';
reconfiles = natdir([recon_folder,'*.tif']);
% flfm_folder = '.\rawtif_selected_r_selected\rawtif_selected_r_1000_selected_offsetsub_acsn_ccut\';
recon_folder = '.\PSFFLFint_20230706_Red_Gly_10um-bead1000lf_ccut_single\';
% recon_folder = '.\bead1000_single\';
reconfiles = [reconfiles;natdir([recon_folder,'*.tif'])];
% flfm_folder = '.\rawtif_selected_r_selected\rawtif_selected_r_2000_selected_offsetsub_acsn_ccut\';
recon_folder = '.\PSFFLFint_20230706_Red_Gly_10um-bead2000lf_ccut_single\';
% recon_folder = '.\bead2000_single\';
reconfiles = [reconfiles;natdir([recon_folder,'*.tif'])];
% flfm_folder = '.\rawtif_selected_r_selected\rawtif_selected_r_4000_selected_offsetsub_acsn_ccut\';
recon_folder = '.\PSFFLFint_20230706_Red_Gly_10um-bead4000lf_ccut_single\';
% recon_folder = '.\bead4000_single\';
reconfiles = [reconfiles;natdir([recon_folder,'*.tif'])];
%%
clc,volstats = {};
h = waitbar(0,'Calculating volumes...');
fig = figure(101);clf,hold on;
fig.NumberTitle = "off";
% for ii = 50%[50,100,200,250]
for ii = 1:1:length(reconfiles)
    reconfilename = [reconfiles(ii).folder,'\',reconfiles(ii).name];
    imrecon = single(tiffreadVolume(reconfilename,'PixelRegion',{[1 inf],[1 inf],[1 101]}));
    imrecon_resize = imresize3(imrecon,[size(imrecon,1) size(imrecon,2) round(size(imrecon,3)*100/65)]);
    imrecon_resize = rescale(imrecon_resize,0,1);
    % imrecon = imrecon.*single(imrecon>0.1);
    imrecon_proj = rescale(mean(single(imrecon),3),0,1);
    
    % -----------------------------------------
    % imrecon_max = rescale(double(squeeze(max(imrecon,[],[1 2]))),0,1);
    
    try
        [Rx,Ry,Rz] = getRxyz(imrecon_resize,3.3,"lf");
        imV = 4/3*pi*Rx*Ry*Rz;
        imD = (imV*6/pi)^(1/3)*0.065;
% ====================================================

%         imrecon_vol = pi/6* prod(Sizes,'all');
        imint = sum(imrecon,"all")/1e8;
        if isempty(Rz)
            continue
        else
            volstats = cat(1,volstats,{reconfiles(ii).name,Rz,imD,imint});
        end
        save('volstats_Red_lf.mat',"volstats")
    catch ME
        % volstats = cat(1,volstats,{reconfiles(ii).name,-1,-1,0});
        % save('volstats_Blue_lf.mat',"volstats")
        disp(['Edge Detect ERR in "',reconfiles(ii).name,'"!'])
        continue
        % rethrow(ME)
    end

    subax1 = subplot(1,2,1,Parent=fig);imshow(imrecon_proj,Parent=subax1);
    subax2 = subplot(1,2,2,Parent=fig);scatter(subax2,cell2mat(volstats(:,end-1)),...
                                                      cell2mat(volstats(:,end)), ...
                                                      14,[0.5 0.5 1],'filled');
    fig.Name = [num2str(ii),' | ',reconfiles(ii).name,' | ',num2str(imD)];
    drawnow;
    waitbar(ii/length(reconfiles),h,[num2str(ii),'/',num2str(length(reconfiles)),': ',reconfiles(ii).name])
end
hold off
delete(h)
disp('Done!')
figure(102),histogram(cell2mat(volstats(:,end-1)),80);
%%
% clc
% for ii = 1:461
%     try
%         if volstats{ii,4}>4.6
%             volstats(ii,:) = [];
%         end
%     catch
%         continue
%     end
% end