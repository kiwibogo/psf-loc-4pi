% for v19 iPALMast_analysisv14scmos.m
% v19 
% includes sCMOS noise

function [rois,cor, T] = makerois_dual(im1,im2,T,thresh,bxsz, updateT)
sz = 16/5;
exROI = [40,0,0,0];
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


% if updateT == 1
%     for kk = 1:1
%         [rois] = getrois(im1,im2,T,cor1,bxsz);
%         roi = squeeze(sum(rois,3));
%         Nchannel = size(roi,3);
%         xi = zeros(1,Nchannel);
%         yi = zeros(1,Nchannel);
%         for ss = 1:Nchannel
%             [yi(ss),xi(ss)] = find(roi(:,:,ss)==max(roi(:,:,ss),[],'all'));
%             
%         end
%         dx = xi-xi(1);
%         dy = yi-yi(1);
%         
%         T(1:2,3,:) = T(1:2,3,:)+reshape([dx;dy],2,1,Nchannel);
%     end
% end

[rois,cor] = getrois(im1,im2,T,cor1,bxsz);
s1 = floor(bxsz/2);
cor(:,1:2,:) = cor(:,1:2,:)-s1; % set coordinate to the upper left corner of the subregion

end

function [rois,cor] = getrois(im1,im2,T,cor1,bxsz)
    rois = [];
    cor = [];
    Num = size(cor1,1);
    imsz = size(im1);
    imgcc = imsz([2,1])/2;
    for ss = 1:2
        eval(['ims = im',num2str(ss),';'])
        %corT = o.transformToTarget(ss,[cor1(:,1),cor1(:,2)]);


        corT = (T(:,:,ss)*[cor1(:,[1,2])-imgcc,ones(Num,1)]')';
        corT = corT(:,1:2)+imgcc;

        [sub1] = cMakeSubregions(round(corT(:,2)),round(corT(:,1)),cor1(:,3),bxsz,single(ims));
        rois = cat(4,rois,sub1);
        cor = cat(3,cor,[corT(:,[1,2]),cor1(:,3)]);
    end

end