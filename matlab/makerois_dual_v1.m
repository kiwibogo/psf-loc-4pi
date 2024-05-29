% for v19 iPALMast_analysisv14scmos.m
% v19 
% includes sCMOS noise

function [rois,cor,roisum] = makerois_dual(im1,im2,T,thresh,bxsz)
sz = 16/5;
exROI = [40,0,0,0];
imsz = size(im1);
corall = cell(1,2);
for ss = 1:2
    eval(['ims = im',num2str(ss),';'])
    ims(ims<0) = min(abs(ims(:)));
    rangemin = ceil(bxsz/2);
    rangemax = imsz(1)-ceil(bxsz/2)+1;
    im_unif=unif(ims,[sz sz 0],'rectangular')-unif(ims,[2*sz 2*sz 0],'rectangular');
    im_max=(im_unif>=.999*maxf(im_unif,[2*sz 2*sz 0],'rectangular'))&(im_unif>thresh);
    centers = findcoord(im_max);
    %     impeak = cHistRecon3D(imsz(1),imsz(2),imsz(3),single(centers(:,2)),single(centers(:,1)),single(centers(:,3)),0);
    %     Imax = max(sumim(:));
    %     overlay(sumim./Imax.*500,impeak)
    cor = centers;
    mask = cor(:,1)<rangemax & cor(:,1)>rangemin & cor(:,2)<rangemax & cor(:,2)>rangemin;
    corall{ss} = cor(mask,:);
    eval(['cor',num2str(ss),'=cor;'])
end
imgcc = imsz([1,2])/2;
corT = (T(:,:,2)\[cor2(:,[1,2])-imgcc,ones(size(cor2,1),1)]')';
corT = corT(:,1:2)+imgcc;

% remove nearby coordinates
xk = cat(2,cor1(:,1),corT(:,1));
yk = cat(2,cor1(:,2),corT(:,2));
fk = cat(2,cor1(:,3),corT(:,3));
dis = pdist([xk,yk]);
link = linkage(dis,'complete');
Tc = cluster(link,'cutoff',bxsz/2,'criterion','distance');

posM = []; % global coordinate
for tt=1:max(Tc)
    maskT=(Tc==tt);
    if sum(maskT)==1
        posM = cat(1,posM,[xk(maskT),yk(maskT),fk(maskT)]);
    end
end


cor1 = posM;
rois = [];
cor = [];
Num = size(cor1,1);
imgcc = imsz([1,2])/2;
for ss = 1:4
    eval(['ims = qd',num2str(ss),';'])
    %corT = o.transformToTarget(ss,[cor1(:,1),cor1(:,2)]);
    

    corT = (T(:,:,ss)*[cor1(:,[1,2])-imgcc,ones(Num,1)]')';
    corT = corT(:,1:2)+imgcc;
    
    [sub1] = cMakeSubregions(round(corT(:,2)),round(corT(:,1)),cor1(:,3),bxsz,single(ims));
    rois = cat(4,rois,sub1);
    cor = cat(3,cor,[corT(:,[2,1]),cor1(:,3)]);
end

s1 = floor(bxsz/2);
cor(:,1:2,:) = cor(:,1:2,:)-s1; % set coordinate to the upper left corner of the subregion

end