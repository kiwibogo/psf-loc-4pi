
function [P,CRLB,LogL] = psfloc(obj,rois,dx,dy,shared)
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
iterations = 100;

P = zeros(numlocs,7+3*(6-sum(shared)));
CRLB = zeros(numlocs,6+3*(6-sum(shared)));
fnum = 10000;
[cormask] = gencormask([0:numlocs-1],[],fnum);
Nf = max(cormask);

dTAll = cat(1,reshape(dx',1,4,[]),reshape(dy',1,4,[]),zeros(1,4,numlocs),zeros(1,4,numlocs),repmat(obj.Dz,1,1,numlocs),repmat(obj.Dphi,1,1,numlocs));

for p=1:sp0(2)
    for k=1:sz0(2)
        Ph = [];CRLBh=[];LogLh=[];
        for kk = 1:Nf
            mask = cormask==kk;
            [ph,ch, Lh, psfh] = GPUmleFit_LM_4Pi_v1(single(data(:,:,mask,:)),uint32(sharedA(:,mask)),iterations,single(obj.IABall),single(dTAll(:,:,mask)),single(phi0A(:,mask)),single(zstart(mask,k)),single(p0(mask,p)),1);
            Ph = cat(1,Ph,ph);
            CRLBh = cat(1,CRLBh,ch);
            LogLh = cat(1,LogLh,Lh);
            clear GPUmleFit_LM_4Pi_v1
        end
        
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end
    
    
end
