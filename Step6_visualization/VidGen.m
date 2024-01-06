clear,clc
addpath("Z:\Xuanwen\FLFMuf\ExpData\utilities\")
%% Set stack sequences of beads
stackfolder = ".\beads\";
stackseq = ["bead_wf_r_00001(3)_20_61",".tif","r","Wide field";
            "bead_wf_g_00001(3)_30_71",".tif","g","Wide field";
            "bead_wf_b_00001(5)_50_91",".tif","b","Wide field"];
%% Set stack sequences of peromito
stackfolder = ".\peromito\";
stackseq = ["pemi_wf00002(1)_240_281",".tif","br","Wide field";
            "pemi_00002(3)_910_951",".tif","br","Light field";
            "V_pemi_00002(3)_910_951",".tif","c","Volume rendering"];
%% Set stack sequences of JurkatSTS
stackfolder = ".\JurkatSTS\";
stackseq = ["JK_minu_wf00001(10)_1_42",".tif","rbbrbbrrbbrbbrb"+"brbrbrbrbrbrbrbrb"+"brbrrrbrbr","Wide field";
            "JK_minu30_wf_00001(3)_6_47",".tif","rbrbbrbrbrbrbrrb"+"rbbrbrbrbrbbrbrrb"+"rbrbrbrbb","Wide field";
            "JK_minu60_wf_00001(10)_43_84",".tif","rbrbrbbrrbrbrbrb"+"brbrbrrbrbrbrb"+"rbbrrbrbrbrb","Wide field";
            "JK_minu120_wf_00001(3)_1_42",".tif","brbrbrbrrbrbrbrb"+"rbrbrbrbrbrbrbrbrbrbrbrbrb","Wide field";
            "JK_minu300_wf_00001(3)_101_142",".tif","rbrbrbrbrbrbrbrbrbrbrbrb"+"brbrbrbrbrbrbrbrbr","Wide field"];
% 
% stackseq = ["JK_minu_00001(2)_71_112",".tif","rbrbrbrbbrbbrbrbbrbrr"+"brbrbbrrbrrbrbr"+"rbrbrb","Light field";
%             "JK_minu30_00001(5)_250_291",".tif","rbrbrbrbrbrb"+"brrbrbrb"+"rbrbrbrbrb"+"brbrbrbrbrbr","Light field";
%             "JK_minu60_00001(2)_145_186",".tif","rbbr"+"brrbrbrbbrbr"+"brrrbrbrbrrbbrb"+"rbrbbbbrbrb","Light field";
%             "JK_minu120_00001(2)_1_42",".tif","brbrbrbrbrbr"+"rbrbbrbrbrbrbr"+"rbrbrbbrbrbrrbrb","Light field";
%             "JK_minu300_00002(1)_8_49",".tif","brbrbrbrbr"+"rbrbrbrbrb"+"brbrrbbrbr"+"brrbrbrbrbrb","Light field"];

% stackseq = ["V_JK_minu_00001(2)_71_112",".tif","c","Volume rendering";
%             "V_JK_minu30_00001(5)_250_291",".tif","c","Volume rendering";
%             "V_JK_minu60_00001(2)_145_186",".tif","c","Volume rendering";
%             "V_JK_minu120_00001(2)_1_42",".tif","c","Volume rendering";
%             "V_JK_minu300_00002(1)_8_49",".tif","c","Volume rendering"];
%% Set stack sequences of VLS
stackfolder = ".\VLS\";
stackseq = ["S6_wf_00001(3)_310_351",".tif","gr","Wide field"];
stackseq = ["L5_wf_00001(2)_870_911",".tif","rg","Wide field"];
%% Set stack sequences of spleen
stackfolder = ".\bloodspleen\";
stackseq = ["spleen_wf_00001(5)_1_42",".tif","b","Wide field";
            "spleen_00001(3)_355_396",".tif","b","Light field";
            "V_spleen_00001(3)_355_396",".tif","c","Volume rendering"];
%% Set stack sequences of blood
stackfolder = ".\bloodspleen\";
stackseq = ["blood_wf_00001(1)_220_261",".tif","b","Wide field";
            "blood_00001(11)_180_221",".tif","b","Light field";
            "V_blood_00001(11)_180_221",".tif","c","Volume rendering"];
%% write wideo
thres = 90;
savingfolder = "C:\Users\Xuanwen\Desktop\videos\";
% savingfolder = ".\vid\";
fontsize = 32;
btmloc = 800;
% =====================================================================
imstack = [];
for ii = 1:size(stackseq,1)
    istack = single(tiffreadVolume(stackfolder+stackseq(ii,1)+stackseq(ii,2)));
    if stackseq(ii,3)~="c"
        istack = (istack-90).*((istack-thres)>=0);
    else
        resizesize = 1024;
        istack = cat(4,imresize3(istack(:,:,:,1),round([resizesize,resizesize,42])),...
                       imresize3(istack(:,:,:,2),round([resizesize,resizesize,42])),...
                       imresize3(istack(:,:,:,3),round([resizesize,resizesize,42])));
    end
    imstack = cat(1,imstack,{istack});
end
%%
clc;
v = VideoWriter(savingfolder+stackseq(1));
v.FrameRate = 10;
open(v)
h = waitbar(0,"Writing videos...");

for ii = 1:size(imstack{1},3)
    ims = [];
    for jj = 1:size(imstack,1)
        istack = imstack{jj};
        colorseq = char(stackseq(jj,3));
        colorii = colorseq(length(colorseq)-mod(length(colorseq)-ii,length(colorseq)));
        if colorii == 'r'
            im = rescale(istack(:,:,ii),0,1);
            im = uint8(cat(3,im,im*0,im)*255);
            im = insertText(im,[50,btmloc],stackseq(jj,4),...
                'Font','Calibri','FontSize',fontsize,'TextColor','white','BoxOpacity',0);
            if stackseq(jj,4)=="Wide field"
                im = drawScalebar(im,1e-5,65e-9,'downright',[0.1,0.02],'Calibri',32,'white');
            else
                im = drawScalebar(im,1e-5,145e-9,'downright',[0.1,0.1],'Calibri',32,'white');
            end
        elseif colorii == 'g'
            im = rescale(istack(:,:,ii),0,1);
            im = uint8(cat(3,im*0,im,im*0)*255);
            im = insertText(im,[50,btmloc],stackseq(jj,4),...
                'Font','Calibri','FontSize',fontsize,'TextColor','white','BoxOpacity',0);
            if stackseq(jj,4)=="Wide field"
                im = drawScalebar(im,1e-5,65e-9,'downright',[0.1,0.02],'Calibri',32,'white');
            else
                im = drawScalebar(im,1e-5,145e-9,'downright',[0.1,0.1],'Calibri',32,'white');
            end
        elseif colorii == 'b'
            im = rescale(istack(:,:,ii),0,1);
            im = uint8(cat(3,im*0,im,im)*255);
            im = insertText(im,[50,btmloc],stackseq(jj,4),...
                'Font','Calibri','FontSize',fontsize,'TextColor','white','BoxOpacity',0);
            if stackseq(jj,4)=="Wide field"
                im = drawScalebar(im,1e-5,65e-9,'downright',[0.1,0.02],'Calibri',32,'white');
            else
                im = drawScalebar(im,1e-5,145e-9,'downright',[0.1,0.1],'Calibri',32,'white');
            end
        elseif colorii == 'c'
            im = squeeze(uint8(istack(:,:,ii,:)));
        end
        ims = cat(2,ims,im,uint8(ones(size(im,1),20,3)*255*(jj~=size(stackseq,1))));
    end
    ims = insertText(ims,[50,50],sprintf('%.3f',round((ii-1)*0.005,3))+" sec",...
            'Font','Calibri','FontSize',fontsize,'TextColor','white','BoxOpacity',0);
    writeVideo(v,ims);
    waitbar(ii/size(imstack,3),h,"Writing videos "+round(ii/size(imstack,3)*100,2)+"% completed!")
end

close(v);
close(h);
