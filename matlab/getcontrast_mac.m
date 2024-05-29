% for v14 iPALMast_analysisv14scmos.m

% Version numbers:
% July30th 2014: now crop each image and their noisemap accordingly in
% their native scmos qd images.

% previous versions could be located from google drive versions.

function [rms, rmp,ang_contrast] = getcontrast_mac(q1,q2,q3,q4,subsz,res)
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
    subims(:,:,:,ii) = pickbead(im,subsz,N,'single',[round(trans_y),round(trans_x),tc],[]);
    xcenter(:,ii) = trans_x(:)-round(trans_x(:));
    ycenter(:,ii) = trans_y(:)-round(trans_y(:));
    
end
sigma = 0.9;
sz = size(subims,1);
N = numel(xf);
models = zeros(sz,sz,N,4);

for ii = 1:4
    out = finitegausspsf(sz,[sigma sigma],1,0,[xcenter(:,ii),ycenter(:,ii)]+floor(sz/2));
    models(:,:,:,ii) = out;
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




