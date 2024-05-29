function [IABallo,F, ccz, shiftxy, T] = getIABmodel(datapath, dataname, Iratio)

F = load([datapath,'\',dataname,'.mat']);
T = double(cat(3,eye(3,3),permute(F.T,[3,2,1])));
shiftxy = squeeze(T([2,1],3,:))';

IL = permute(F.I_model,[3,4,2,1]);
AL = permute(F.A_model,[3,4,2,1]);
IABallo = [];
if nargin<3
    normf = squeeze(sum(sum(IL(:,:,:,1)-F.offset,1),2));
    normf = max(normf)*2*ones(1,4);
else
    normf = squeeze(sum(sum(IL(:,:,:,1)-F.offset,1),2));
    normf = max(normf)*2./Iratio;
end
for ss = 1:4
    I = IL(:,:,:,ss);
    A = 2*real(AL(:,:,:,ss));
    B = -2*imag(AL(:,:,:,ss));
    
    Ism = I; Asm = A; Bsm = B;
%     lambdax = 0.00;
%     lambdaz = 0.2;
%     lambda = [lambdax lambdax lambdaz];
%     Ism = getsmoothpsf(I,lambda);
%     Asm = getsmoothpsf(A,lambda);
%     Bsm = getsmoothpsf(B,lambda);
    
    %imcom = cat(3,Ism,Ism-Asm,Ism+Asm,Ism-Bsm,Ism+Bsm);
    %offset = min(imcom(:));
    offset = F.offset;
    I1 = Ism-offset;
    
    sz = size(I1);
    ccz = floor(sz(3)/2);
    cc = ceil(sz(1)/2);
    bx = cc-1;
    
    I2 = I1./normf(ss);
    A2 = Asm./normf(ss);
    B2 = Bsm./normf(ss);
    
    Ispline = Spline3D_interp(I2);
    Aspline = Spline3D_interp(A2);
    Bspline = Spline3D_interp(B2);
    
    IAB = cat(5,Aspline,Bspline,Ispline);
    
    IABallo = cat(6,IABallo,IAB);
end

indx = cc+1; indy = cc+1;
h = figure('position',[200,200,950,750]);
subplot(221)
plot([squeeze(abs(I(indx,indy,:))),squeeze(abs(Ism(indx,indy,:)))])
title('I model')
subplot(222)
plot([squeeze((A(indx,indy,:))),squeeze((Asm(indx,indy,:)))])
title('A model')
subplot(223)
plot([squeeze(abs(I2(indx,indy,:))),squeeze(abs(Ism(indx,indy,:)))./normf])
title('I model')
subplot(224)
plot([squeeze((A2(indx,indy,:))),squeeze((Asm(indx,indy,:)))./normf])
title('A model')


h = figure('position',[200,300,1100,300]);
subplot(131)
plot(squeeze(sum(IL,[1,2]))./normf,'linewidth',1.2); axis tight
subplot(132)
plot(squeeze(sum(2*abs(AL),[1,2]))./normf,'linewidth',1.2); axis tight
subplot(133);hold on
g = 2*abs(AL)./IL;
plot(squeeze(g(indx,indy,:,:)),'linewidth',1.2)

Isum = squeeze(sum(IL,[1,2]))./normf;
Asum = squeeze(sum(2*abs(AL),[1,2]))./normf;

Isp = (Isum(:,1)+Isum(:,3))./(Isum(:,2)+Isum(:,4));
Asp = (Asum(:,1)+Asum(:,3))./(Asum(:,2)+Asum(:,4));
figure;plot([Isp,Asp]);
legend('p/s I','p/s A')


