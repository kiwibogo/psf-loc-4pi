function [PSFs_blur,zemit] = vectorPSF_single(Fm,params, zrange, stagepos, fillfactor, shift, psfcenter)
% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), sampling in pupil
%optics para
% add norm function comparecd to v3 #need test
% add pixelX Y


NA = params.option.imaging.NA;
refmed = params.option.imaging.RI.med;
refcov = params.option.imaging.RI.cov;
refimm = params.option.imaging.RI.imm;
lambda = params.option.imaging.emission_wavelength;

z_center = stagepos*refmed/refimm;
if zrange(1)<-z_center
    zrange = zrange-zrange(1)-z_center;
end
%image para
%zemit0 = parameters.zemit0;
xemit = zeros(size(zrange));
yemit = zeros(size(zrange));
zemit = zrange+z_center;
pupil = permute(Fm.res.pupil,[1,2]);
Npupil = size(pupil,1);
I_model = permute(Fm.res.I_model,[3,2,1]);
sizeX = params.roi.roi_size(1);
sizeY = params.roi.roi_size(2);
sizeZ = numel(zemit);
pixelSizeX = params.pixel_size.x;
pixelSizeY = params.pixel_size.y;
xrange = pixelSizeX*sizeX/2;
yrange = pixelSizeY*sizeY/2;

% pupil radius (in diffraction units) and pupil coordinate sampling
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid
CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
kx = (2*pi*NA/lambda)*XPupil;
ky = (2*pi*NA/lambda)*YPupil;
kz = (2*pi*refimm/lambda)*CosThetaImm;
kzmed = (2*pi*refmed/lambda)*CosThetaMed;

%need check again, compared to original equation, FresnelPmedcov is
%multipiled by refmed
FresnelPmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% Apoidization for sine condition
%apoid = sqrt(CosThetaImm)./CosThetaMed;
apoid = 1;
% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
if nargin>4
     
    ApertureU = double((XPupil.^2+YPupil.^2)<fillfactor);
    shiftphase = exp(1i*(shift*pixelSizeX)*kx.*ApertureU);
    pupil = pupil.*shiftphase;
    switch psfcenter
        case 'ua'
           xemit = xemit+shift*pixelSizeX;
        case 'sa'
            xemit = xemit;
    end

Amplitude = ApertureMask.*apoid;

% setting of vectorial functions
Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

pvec{1} = FresnelP.*CosTheta.*CosPhi;
pvec{2} = FresnelP.*CosTheta.*SinPhi;
pvec{3} = -FresnelP.*SinTheta;
svec{1} = -FresnelS.*SinPhi;
svec{2} = FresnelS.*CosPhi;
svec{3} = 0;

PolarizationVector = cell(2,3);
for jtel = 1:3
    PolarizationVector{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
    PolarizationVector{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
end



% calculation aberration function
% if isempty(parameters.pupil)
%     Waberration = zeros(size(XPupil));
%     orders = parameters.aberrations(:,1:2);
%     zernikecoefs = squeeze(parameters.aberrations(:,3));
%     normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
%     zernikecoefs = normfac.*zernikecoefs;
%     allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
%     for j = 1:numel(zernikecoefs)
%       Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));
%     end
%
%     pupil = exp(2*pi*1i*Waberration/lambda);
% else
%     pupil = parameters.pupil;
% end

PupilMatrix = cell(2,3);
for itel = 1:2
    for jtel = 1:3
        PupilMatrix{itel,jtel} = Amplitude.*pupil.*PolarizationVector{itel,jtel};
    end
end

% pupil and image size (in diffraction units)
% PupilSize = NA/lambda;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

% calculate auxiliary vectors for chirpz
[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,sizeX);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,sizeY);

FieldMatrix = cell(2,3,sizeZ);
for jz = 1:sizeZ
    phixy = -xemit(jz)*kx-yemit(jz)*ky;
    phiz = kzmed*zemit(jz)-kz*stagepos;
    % xyz induced phase contribution

    PositionPhaseMask = exp(1i*(phixy+phiz));


    for itel = 1:2
        for jtel = 1:3

            % pupil functions and FT to matrix elements
            PupilFunction =PositionPhaseMask.*PupilMatrix{itel,jtel};
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
        end
    end
end

%calculates the free dipole PSF given the field matrix.
PSFs = zeros(sizeX,sizeY,sizeZ);
for jz = 1:sizeZ
    for jtel = 1:3
        for itel = 1:2
            PSFs(:,:,jz) = PSFs(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
        end
    end
end

% calculate intensity normalization function using the PSFs at focus
% position without any aberration. It might not work when sizex and sizeY
% are very small
%FieldMatrix = cell(2,3);
pupilsum = 0;
for itel = 1:2
    for jtel = 1:3
        PupilFunction = Amplitude.*PolarizationVector{itel,jtel};
        pupilsum = pupilsum+ sum(abs(PupilFunction).^2,"all");
        %IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
        %FieldMatrix{itel,jtel} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
    end
end

% intFocus = zeros(sizeX,sizeY);
% for jtel = 1:3
%     for itel = 1:2
%         intFocus = intFocus + (1/3)*abs(FieldMatrix{itel,jtel}).^2;
%     end
% end

normIntensity = sqrt(pupilsum/6);

PSFs = PSFs./normIntensity;

% h = fspecial('gaussian',5,0.5);
% PSFs = convn(PSFs,h,'same');

sigma = Fm.res.sigma*pi;
xrange = linspace(-sizeX/2+0.5,sizeX/2-0.5,sizeX);
[xx,yy] = meshgrid(xrange,xrange);
pkx = xx/sizeX;
pky = yy/sizeX  ;
kspace_x = pkx.*pkx;
kspace_y = pky.*pky;

filter2 = exp(-2*sigma(2)*sigma(2)*kspace_x-2*sigma(1)*sigma(1)*kspace_y);

filter2 = filter2/max(filter2);
PSFs_blur = zeros(sizeX,sizeY,sizeZ);
for jz = 1:sizeZ
    otf = fftshift(fft2(fftshift(PSFs(:,:,jz))));
    PSFs_blur(:,:,jz) = ifftshift(ifft2(ifftshift(otf.*filter2)));

end


end





