clear,clc
foldname = './PSF_FLFM100X_downsample/';
prefix = 'PSF_';

for fileind = -5.0:0.1:5.0
    load([foldname,prefix,num2str(fileind,'%1.1f'),'.mat'],'PSF')
    PSF = PSF(1:2:6001,1:2:6001);
    save([foldname,prefix,num2str(fileind,'%1.1f'),'.mat'],'PSF')
    disp([foldname,prefix,num2str(fileind,'%1.1f'),'.mat'])
end