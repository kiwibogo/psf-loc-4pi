function [P,CRLB,LogL,PSF] = psfloc_ast_single(obj,data)
numlocs = size(data,3);

zstart = obj.Initz;

varmap = 0;
silent = 0;

LogL = ones(numlocs,1,'single')-1e10;

P = zeros(numlocs,6);
CRLB = zeros(numlocs,5);
fnum = 20000;
[cormask] = gencormask([0:numlocs-1],[],fnum);
Nf = max(cormask);

fittype = 5;
iterations = 100;

datasize = size(data);
PSF = zeros(datasize);


for k=1:size(zstart,2)
    Ph = [];CRLBh=[];LogLh=[];PSFh = [];
    for kk = 1:Nf
        mask = cormask==kk;
        [ph,ch, Lh,psfh]=GPUmleFit_LM_v1(single(data(:,:,mask)),fittype,iterations,single(obj.coeff),varmap,silent,single(zstart(k)));
        
        %[ph,ch, Lh]=GPUmleFit_LM(single(data(:,:,mask)),fittype,iterations,single(obj.coeff),varmap,silent,single(zstart(k)));

        datsz = size(data(:,:,mask));
        Ph = cat(1,Ph,ph);
        CRLBh = cat(1,CRLBh,ch);
        LogLh = cat(1,LogLh,Lh);
        PSFh = cat(3,PSFh,reshape(psfh,datsz));
        clear GPUmleFit_LM_v1
    end
    indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
    P(indbetter,:)=Ph(indbetter,:);
    CRLB(indbetter,:)=CRLBh(indbetter,:);
    LogL(indbetter)=LogLh(indbetter);
    PSF(:,:,indbetter,:) = PSFh(:,:,indbetter,:);

end

