function selpath = CutCircleSave(selpath,boundmat,circlemat)
%CUTSAVE Summary of this function goes here
%   Detailed explanation goes here
filepath = dir([selpath,'\*.tif']);
filenames = string({filepath(:).name}');
mkdir([selpath,'_ccut\'])
for i = 1:length(filenames)
    try
        img = imread([selpath,'\',char(filenames(i))]);
        img = uint16(double(img) .* double(circlemat)+3);
        imgout = img(boundmat(1):boundmat(2),boundmat(3):boundmat(4));
        imwrite(imgout,[selpath,'_ccut\',char(filenames(i))])
        disp(['Saved! (.\..._ccut\',char(filenames(i))])
    catch ME
        disp(ME)
        disp(['Image Unsaved :-( ',num2str(i)])
    end
end

end

