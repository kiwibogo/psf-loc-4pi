% for v19 iPALMast_analysisv14scmos.m
% v19 
% includes sCMOS noise

function [rois,cor1] = makerois(im1,thresh,bxsz)
sz = 5;
exROI = [0,0,0,0];
imsz = size(im1);

ims = im1;
ims(ims<0) = min(abs(ims(:)));
rangemin = ceil(bxsz/2);
im_unif=unif(ims,[sz sz 0],'rectangular')-unif(ims,[2*sz 2*sz 0],'rectangular');
im_max=(im_unif>=.999*maxf(im_unif,[2*sz 2*sz 0],'rectangular'))&(im_unif>thresh);
centers = findcoord(im_max);
%     impeak = cHistRecon3D(imsz(1),imsz(2),imsz(3),single(centers(:,2)),single(centers(:,1)),single(centers(:,3)),0);
%     Imax = max(sumim(:));
%     overlay(sumim./Imax.*500,impeak)
cor = centers;
mask = cor(:,1)<imsz(2)-rangemin+1-exROI(4) & cor(:,1)>rangemin+exROI(3) & cor(:,2)<imsz(1)-rangemin+1-exROI(2) & cor(:,2)>rangemin+exROI(1);
cor1 = cor(mask,:);


[rois] = cMakeSubregions(round(cor1(:,2)),round(cor1(:,1)),cor1(:,3),bxsz,single(ims));


s1 = floor(bxsz/2);
cor1(:,1:2) = cor1(:,1:2)-s1; % set coordinate to the upper left corner of the subregion

end

