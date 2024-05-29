% for v19 iPALMast_analysisv14scmos.m
% v19 
% includes sCMOS noise

function [rois,cor,roisum] = makerois_4pi_mac(qd1,qd2,qd3,qd4,T,thresh,bxsz,shiftxy)
sz = 16/5;

imsz = size(qd1);
sumim = qd1+imtranslate(qd3,shiftxy);

sumim(sumim<0) = min(abs(sumim(:)));
rangemin = bxsz;
rangemax = imsz(1)-bxsz;
im_unif=unif(sumim,[sz sz 0],'rectangular')-unif(sumim,[2*sz 2*sz 0],'rectangular');
im_max=(im_unif>=.999*maxf(im_unif,[2*sz 2*sz 0],'rectangular'))&(im_unif>thresh);
centers = findcoord(im_max);
cor1 = centers;
mask = cor1(:,1)<rangemax & cor1(:,1)>rangemin & cor1(:,2)<rangemax & cor1(:,2)>rangemin;
cor1 = cor1(mask,:);

N = size(cor1,1);
[roisum] = pickbead(sumim,bxsz,N,'single',cor1,[]);

rois = [];
cor = [];
Num = size(cor1,1);
imgcc = imsz([1,2])/2;
for ss = 1:4
    eval(['ims = qd',num2str(ss),';'])

    corT = (T(:,:,ss)*[cor1(:,[1,2])-imgcc,ones(Num,1)]')';
    corT = corT(:,1:2)+imgcc;
    [sub1] = pickbead(ims,bxsz,N,'single',round([corT(:,[1,2]),cor1(:,3)]),[]);

    rois = cat(4,rois,sub1);
    cor = cat(3,cor,[corT(:,[2,1]),cor1(:,3)]);
end

s1 = floor(bxsz/2);
cor(:,1:2,:) = cor(:,1:2,:)-s1; % set coordinate to the upper left corner of the subregion

end