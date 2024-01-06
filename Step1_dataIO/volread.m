function imstack = volread(stack1_name)
%VOLREAD Summary of this function goes here
%   Detailed explanation goes here
% stack1_name = [vfolder,'pero200hz(33)1iter50.tif'];
ii = 1;
imstack = [];
disp('Reading...')
while 1
%     fprintf(['\b|',num2str(ii),'\n']);
    try
        im = imread(stack1_name,'index',ii);
        imstack = cat(3, imstack,im);
        ii = ii + 1;
    catch ME
        if size(imstack,3)<1
            rethrow(ME)
        else
%             disp(stack1_name);
%             disp(['The stack size is [',num2str(size(imstack)),'].']);
            break
        end
    end
end
end

