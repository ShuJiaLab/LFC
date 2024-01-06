function dcimg2tiff(scrnamelist,destfolder,varargin)
%DCIMG2TIFF Summary of this function goes here
%   Detailed explanation goes here

mod = py.importlib.import_module('dcimg');
py.importlib.reload(mod);
if exist(destfolder,'dir')==0
    mkdir(destfolder)
end

for scrind = 1:1:size(scrnamelist,1)
    [~,name,ext] = fileparts(scrnamelist(scrind));
    destname = [destfolder,'\',char(name)];
    if ext ~= ".dcimg"
        disp(['[ERR] File \"',char(name),char(ext),'\" should be in *.dcimg format.'])
        break
    else
        disp(['[NOTE] Files will be saved in ',destname,'(*).tif']);
        scrfile = py.dcimg.DCIMGFile(py.str(scrnamelist(scrind)));
        scrfile_shape = getProperty(scrfile,"shape");
        tifnum = round(ceil(2^((log2(getProperty(scrfile,"zsize"))+...
                                log2(getProperty(scrfile,"ysize"))+...
                                log2(getProperty(scrfile,"xsize")))-31)));
		pagenum = round(floor(2^(31-log2(getProperty(scrfile,"ysize"))-...
                                    log2(getProperty(scrfile,"xsize")))))-1;
        disp(['File in shape (',num2str(scrfile_shape),'), tif numbers: ',...
            num2str(tifnum),', page numbers: ',num2str(pagenum)])

        for tifind = 1:1:tifnum
            destname_ind = [destname, '(',num2str(tifind),').tif'];
            for pageind  = 1:1:pagenum
                if pageind ==1
                    t = Tiff(destname_ind,'w');
                else
                    t = Tiff(destname_ind,'a');
                end
                t.setTag('ImageLength', scrfile_shape(1));
                t.setTag('ImageWidth',  scrfile_shape(2));
                t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
                t.setTag('Compression', Tiff.Compression.None);
                t.setTag('BitsPerSample', 16);
                t.setTag('SamplesPerPixel', 1);
                t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
                t.write( uint16(pyrun(['scrdata = pyscrfile[',num2str(pagenum*(tifind-1)+pageind),',:,:]'],...
                    "scrdata",pyscrfile=scrfile)) );
                t.close;
                disp([scrind,tifind,pageind])
            end
        end

    end
end
end

function propV = getProperty(dcimgfile,prop)
DCIMG_PROP = ["byte_depth","double";
              "bytes_per_img","double";
              "bytes_per_row","double";
              "deep_copy_enabled","";
              "dtype","";
              "file_size","double";
              "framestamps","double";
              "nfrms","double";
              "shape","cell";
              "timestamps","double";
              "xsize","double";
              "ysize","double";
              "zsize","double";
              "y0","double";
              "file_path","string";
              "first_4px_correction_enabled","";
              "binning","";
              "x0","double";
              "mma","double";
              "fmt_version","double";
              "mm",""];
if ismember(string(prop),DCIMG_PROP)
    [proprow,~] = find(DCIMG_PROP==string(prop));
    eval(['pypropV = dcimgfile.',char(prop),';'])
    eval(['propV = ',char(DCIMG_PROP(proprow,2)),'(pypropV);'])
    if prop=="shape"
        propV = [double(propV{2}),double(propV{3}),double(propV{1})];
    end
else
    disp([char(prop),' is not accessible!'])
end
end