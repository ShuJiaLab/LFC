clear,clc
addpath("Z:\Xuanwen\FLFMuf\ExpData\utilities\bfmatlab\")
for ii = [4:8,10]
cellidx = ii;
dataset_name = ".\cell"+cellidx+"\cell"+cellidx+"_Structured Illumination.czi";

a = bfopen(char(dataset_name)); % open datafile
b = a{1};                                                  % extract variable with image data
allimages_in = cell2mat(permute(b(:,1),[3 2 1]));          % extract image data
clearvars a b
%
allimages_mito_SR = allimages_in(:,:,1:4:end);
allimages_mito_WF = allimages_in(:,:,2:4:end);
allimages_pero_SR = allimages_in(:,:,3:4:end);
allimages_pero_WF = allimages_in(:,:,4:4:end);
%
output = ".\cell"+cellidx+"\SIM\";
mkdir(output);
output_mitoSR = output + "cell"+cellidx + "_SIM_mito_SR.tif";
output_mitoWF = output + "cell"+cellidx + "_SIM_mito_WF.tif";
output_peroSR = output + "cell"+cellidx + "_SIM_pero_SR.tif";
output_peroWF = output + "cell"+cellidx + "_SIM_pero_WF.tif";
volwrite(allimages_mito_SR,output_mitoSR);
volwrite(allimages_mito_WF,output_mitoWF);
volwrite(allimages_pero_SR,output_peroSR);
volwrite(allimages_pero_WF,output_peroWF);
disp("Saved "+output);
end