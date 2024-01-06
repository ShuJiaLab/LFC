function imV = calcVol_mb(imstack,startP,endP)
%CALVOL_MB Summary of this function goes here
%   Detailed explanation goes here

imV = 0;
H = waitbar(0,'Calculating volume...');
for depth = startP:1:endP
    imlayer = imstack(:,:,depth);
    [~,threshold] = edge(imlayer,'sobel');
    fudgeFactor = 2;
    BWs = edge(imlayer,'sobel',threshold * fudgeFactor);
    strelsize = 5;
    se135 = strel('line',strelsize,135);
    se90 = strel('line',strelsize,90);
    se45 = strel('line',strelsize,45);
    se0 = strel('line',strelsize,0);
    BWsdil = imdilate(BWs,[se135 se90 se45 se0]);
    BWdfill = imclose(BWsdil,strel('disk',30));
    BWnobord = imclearborder(BWdfill,4);
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    imV = imV + sum(BWfinal,'all');                                               
    waitbar((depth-startP)/(endP-startP),H,['Calculating ',num2str(depth),...
        '/',num2str(endP),', (',num2str(round((depth-startP)/(endP-startP),2)),'%)'])
end
delete(H)
end

