function sr = resetparam_imm(params,bxsz,coeff)
pz = params.pixel_size.z;
pxsz = params.pixel_size.x;
Nzm = size(coeff,3)+1;
ccz = size(coeff,3)/2;
sr.Peakthresh = 50;

sr.Boxsize = bxsz;
sr.Gain = 0.5;
sr.ccdoffset = 100;
%sr.Initz = [-1,1]*150/pz+ccz;
sr.Initz = linspace(-Nzm*pz/2,Nzm*pz/2,floor(Nzm*pz/0.5))*0.8/pz+ccz;
sr.Initx = [];
sr.Dz = [0,0];
sr.Iratio = [1,1];
sr.Zcenter = ccz;
sr.coeff = coeff;
sr.Pixelsizex = pxsz*1e3;%nm
sr.Pixelsizey = pxsz*1e3;
sr.Pixelsizez = pz*1e3; % nm

