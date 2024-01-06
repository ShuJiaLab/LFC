clear,clc
% ================ testing ===================
filename = '..\XW20220226\actinTracker_staintest00001.dcimg';

% filename = '\\143.215.226.127\e\XW20220320\bead_flf_r00007.dcimg';
dcimg_hdr_format = {'uint8',[8 1],'file_format';...
                    'uint32',[1 1],'format_version';...
                    'uint32',[5 1],'skip';...
                    'uint32',[1 1],'nsess';...
                    'uint32',[1 1],'nfrms';...
                    'uint32',[1 1],'header_size';...
                    'uint32',[1 1],'skip2';...
                    'uint64',[1 1],'file_size';...
                    'uint32',[2 1],'skip3';...
                    'uint64',[1 1],'file_size2'};
a = memmapfile(filename,'Format',dcimg_hdr_format,'Repeat',1);

dcimg_hdr = a.Data;
dcimg_hdr.file_format = char(dcimg_hdr.file_format');
disp(dcimg_hdr)
% read session header
dcimg_sess_hdr_format = {'uint64',[1 1],'session_size';...
                         'uint32',[13 1],'skip1';...
                         'uint32',[1 1],'nfrms';...
                         'uint32',[1 1],'byte_depth';...
                         'uint32',[1 1],'skip2';...
                         'uint32',[1 1],'xsize';...
                         'uint32',[1 1],'ysize';...
                         'uint32',[1 1],'bytes_per_row';...
                         'uint32',[1 1],'bytes_per_img';...
                         'uint32',[2 1],'skip3';...
                         'uint32',[1 1],'offset_to_data'};
b = memmapfile(filename,'Format',dcimg_sess_hdr_format,'Offset',dcimg_hdr.header_size,'Repeat',1);
dcimg_sess_hdr = b.Data;
disp(dcimg_sess_hdr)
%
dcimg_crop_info_format =  {'uint16',[1 1],'x0';...
                           'uint16',[1 1],'xsize';...
                           'uint16',[1 1],'y0';...
                           'uint16',[1 1],'ysize'};
ii = dcimg_hdr.header_size + 712;
c = memmapfile(filename,'Format',dcimg_crop_info_format,'Offset',ii,'Repeat',1);
dcimg_crop_info = c.Data;
disp(dcimg_crop_info);

%%
clc;
bd = dcimg_sess_hdr.byte_depth;
data_offset = dcimg_hdr.header_size + dcimg_sess_hdr.offset_to_data;
strides = round(double([bd,bd,(dcimg_sess_hdr.bytes_per_img + 32)])/2);
dc_4px_shape = double([1,4,dcimg_sess_hdr.nfrms]);
dc_4px_format = {'uint16',[1 dc_4px_shape(end)*strides(end)],'dc_4px'};
% dc_4px_mem = memmapfile(filename,'Format',dc_4px_format,...
%     'Offset',data_offset + dcimg_sess_hdr.bytes_per_img + 12,'Repeat',1);
dc_4px_mem = memmapfile(filename,'Format',dc_4px_format,...
    'Offset',data_offset,'Repeat',1);
slices = [990:1000];
for ii = slices
    dc_4px_ind = repmat((1:dc_4px_shape(2))',1,dc_4px_shape(1),length(slices))+repmat(0:dc_4px_shape(1)-1,dc_4px_shape(2),1,length(slices))*strides(2);
    dc_4px_ind = dc_4px_ind + repmat(reshape(slices-1,1,1,[]),dc_4px_shape(2),dc_4px_shape(1),1)*strides(3);
end
dc_4px_ind = dc_4px_ind + double(dcimg_sess_hdr.bytes_per_img + 12)/2;
dc_4px = dc_4px_mem.Data.dc_4px(dc_4px_ind(:));
dc_4px = reshape(dc_4px,dc_4px_shape(2),dc_4px_shape(1),length(slices));
dc_4px =permute(dc_4px,[2 1 3]);
disp(squeeze(dc_4px_ind(:,:,1:10))')
disp(squeeze(dc_4px)')
%%
clc;
padding = dcimg_sess_hdr.bytes_per_img - dcimg_sess_hdr.xsize * dcimg_sess_hdr.ysize * bd;
padding = floor(padding/dcimg_sess_hdr.ysize);
data_strides = double(round([bd,dcimg_sess_hdr.xsize*bd+padding,dcimg_sess_hdr.bytes_per_img+32]/2));
data_shape = double([dcimg_sess_hdr.ysize,dcimg_sess_hdr.xsize,dcimg_sess_hdr.nfrms]);
data_size = [1, data_shape(end)*data_strides(end)];
data_format = {'uint16',data_size,'data'};
mma_mem = memmapfile(filename,'Format',data_format,'Offset',data_offset,'Repeat',1);
% slices = 1000;
for ii = slices
    mma_ind = repmat((1:data_shape(2))',1,data_shape(1),length(slices)) + repmat(0:data_shape(1)-1,data_shape(2),1,length(slices))*data_strides(2);
    mma_ind = mma_ind + repmat(reshape(slices-1,1,1,[]),data_shape(2),data_shape(1),1)*data_strides(3);
end
disp('Done!')
%%
clc,
mma = mma_mem.Data.data(mma_ind(:));
mma = reshape(mma,[data_shape(2),data_shape(1),length(slices)]);
mma =permute(mma,[2 1 3]);
disp(mma(1023 - round(dcimg_crop_info.y0)+1,1:4,1))
imshow(double(mma(:,:,1))/max(max(double(mma(:,:,1))))*50)