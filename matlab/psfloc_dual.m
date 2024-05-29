
function [P,CRLB,LogL] = psfloc_dual(obj,rois,dx,dy,shared)
data = rois;
numlocs = size(rois,3);
Nchannel = size(rois,4);
Nparam = 5;

zstart = ones(numlocs,1,'single')*[obj.Initz];
sharedA = repmat(shared',[1 numlocs]);

LogL = ones(numlocs,1,'single')-1e10;
sz0 = size(zstart);
iterations = 100;
param_shift = zeros(Nparam,Nchannel);
param_ratio = ones(Nparam,Nchannel);
param_shift(3,:) = obj.Dz;
param_ratio(4,:) = obj.Iratio;

P = zeros(numlocs,(Nparam+1)+(Nchannel-1)*(Nparam-sum(shared)));
CRLB = zeros(numlocs,Nparam+(Nchannel-1)*(Nparam-sum(shared)));
fnum = 10000;
[cormask] = gencormask([0:numlocs-1],[],fnum);
Nf = max(cormask);
fittype =2;

dTS = zeros(Nparam,Nchannel*2,numlocs);
dTS(1,1:Nchannel,:) = reshape(dx',1,Nchannel,[]);
dTS(2,1:Nchannel,:) = reshape(dy',1,Nchannel,[]);
dTS(:,1:Nchannel,:) = dTS(:,1:Nchannel,:) + reshape(param_shift,Nparam,Nchannel,1);
dTS(:,Nchannel+1:end,:) = repmat(param_ratio,1,1,numlocs);
varmap = single([0]);
silent = single(1);

for k=1:sz0(2)
    Ph = [];CRLBh=[];LogLh=[];
    for kk = 1:Nf
        mask = cormask==kk;
        [ph,ch, Lh] = GPUmleFit_LM_MultiChannel(single(data(:,:,mask,:)),fittype,uint32(sharedA(:,mask)),iterations,single(obj.coeff),single(dTS(:,:,mask)),varmap,silent,single(zstart(mask,k)));
        Ph = cat(1,Ph,ph);
        CRLBh = cat(1,CRLBh,ch);
        LogLh = cat(1,LogLh,Lh);
        clear GPUmleFit_LM_MultiChannel
    end
    
    indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
    P(indbetter,:)=Ph(indbetter,:);
    CRLB(indbetter,:)=CRLBh(indbetter,:);
    LogL(indbetter)=LogLh(indbetter);
end
    
    
end
