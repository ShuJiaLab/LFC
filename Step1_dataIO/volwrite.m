function volwrite(Xguess_norm_re,reconfilename)
Depth_Size = size(Xguess_norm_re,3);
t = Tiff(reconfilename,'w');
tagstruct.ImageLength = size(Xguess_norm_re,1);
tagstruct.ImageWidth = size(Xguess_norm_re,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
h = waitbar(0,'Writing tiff file...');
for depth_index  = 1:Depth_Size
    setTag(t,tagstruct);
    write(t,uint16(Xguess_norm_re(:,:,depth_index)));
    writeDirectory(t);
    waitbar(depth_index/Depth_Size,h,['Writing tiff file... (',num2str(round(depth_index/Depth_Size*100,2)),'%)'])
end
close(t);
delete(h);
end