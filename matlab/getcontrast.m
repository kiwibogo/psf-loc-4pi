% for v14 iPALMast_analysisv14scmos.m

% Version numbers:
% July30th 2014: now crop each image and their noisemap accordingly in
% their native scmos qd images.

% previous versions could be located from google drive versions.

function [rms, rmp,ang_contrast] = getcontrast(q1,q2,q3,q4,subsz,res,dx,dy)
xf = res.x;
yf = res.y;
cor = res.cor;
bg_f = res.bg;

tc = double(cor(:,3));
N = numel(tc);
subims = zeros(subsz,subsz,N,4);
xcenter = zeros(N,4);
ycenter = zeros(N,4);

for ii=1:4
    eval(['im=q' num2str(ii) ';']);
    trans_x = double(cor(:,1,ii)+xf);
    trans_y = double(cor(:,2,ii)+yf);
    imsz = size(im);
    %impeak = cHistRecon3D(imsz(1),imsz(2),imsz(3),single(trans_x(:)),single(trans_y(:)),single(tc),0);
    %overlay(im./max(im(:)).*500,impeak)
    [subims(:,:,:,ii), corx, cory] = cMakeSubregions(round(trans_x(:)),round(trans_y(:)),tc,subsz,single(im));
    xcenter(:,ii) = trans_x(:)-corx;
    ycenter(:,ii) = trans_y(:)-cory;
    
end
sigma = 0.9;
sz = size(subims,1);
N = numel(xf);
models = zeros(sz,sz,N,4);
maxN = 5000;
repN = floor(N/maxN);
vec = cat(1,reshape(repmat([1:repN],maxN,1),repN*maxN,1),(repN+1).*ones(N-repN*maxN,1));
for ii = 1:4
    model = [];
    for nn = 1:repN+1
        maski = vec==nn;
        sampN = sum(maski);
        modeli=GPUgenerateBlobs(single(sz*sampN),single(sz),single(xcenter(maski,ii)+[0:sampN-1]'.*sz),single(ycenter(maski,ii)),single(ones(sampN,1)),single(ones(sampN,1).*sigma),single(ones(sampN,1).*sigma),single(zeros(sampN,1)),1);
        tmp = reshape(modeli',sz,sz,sampN);
        tmp = permute(tmp,[2,1,3]);
        model = cat(3,model,tmp);
    end
models(:,:,:,ii) = model;
end

bgs = reshape(repmat(bg_f',sz*sz,1),sz,sz,N);
int_est = zeros(N,4);

for ii = 1:4
    nom=(subims(:,:,:,ii)-bgs./2).*models(:,:,:,ii);
    denom=models(:,:,:,ii).^2;
    int_est(:,ii)=squeeze(sum(sum(nom,1),2))./squeeze(sum(sum(denom,1),2))./2;% devided by 2 is possibly not necessary
end

rms=(int_est(:,1)-int_est(:,3))./(int_est(:,1)+int_est(:,3));% normalize to the photon count of s polarization
rmp=(int_est(:,4)-int_est(:,2))./(int_est(:,4)+int_est(:,2));% normalize to the photon count of p polarization


ang_contrast=sqrt(rms.^2+rmp.^2);




