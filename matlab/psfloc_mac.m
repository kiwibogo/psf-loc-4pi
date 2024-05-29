
function [P,CRLB,LogL,PSF] = psfloc_mac(obj,rois,dx,dy,shared)
data = rois;
data = data+obj.BGoffset;
numlocs = size(rois,3);

zstart = ones(numlocs,1,'single')*[obj.Initz];
p0 = ones(numlocs,1,'single')*single(obj.InitPhase);
sharedA = repmat(shared',[1 numlocs]);
phi0 = obj.Phi0;
phi0A = repmat(phi0',[1 numlocs]);

LogL = ones(numlocs,1,'single')-1e10;
sz0 = size(zstart);
sp0 = size(p0);
iterations = 10;

P = zeros(numlocs,7+3*(6-sum(shared)));
CRLB = zeros(numlocs,6+3*(6-sum(shared)));

dTAll = cat(1,reshape(dx',1,4,[]),reshape(dy',1,4,[]),zeros(1,4,numlocs),zeros(1,4,numlocs),repmat(obj.Dz,1,1,numlocs),repmat(obj.Dphi,1,1,numlocs));
silent = int32(1);
datasize = size(data);
PSF = zeros(datasize);
for p=1:sp0(2)
    for k=1:sz0(2)
        [Ph,CRLBh, LogLh,psfh] = CPUmleFit_LM_4Pi_v1(single(data),uint32(sharedA),iterations,single(obj.IABall),single(dTAll),single(phi0A),single(zstart(:,k)'),single(p0(:,p)'),silent);
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
        psf = reshape(psfh,datasize);
        PSF(:,:,indbetter,:) = psf(:,:,indbetter,:);
    end
end



clear CPUmleFit_LM_4Pi_v1 
  
    

