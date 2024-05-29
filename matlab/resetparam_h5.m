function sr = resetparam_h5(offsetfile,Fm,params,bxsz,trans)
coeff = permute(Fm.locres.coeff,[6,5,4,3,2,1]);
Tm = double(cat(3,eye(3,3),permute(Fm.res.T,[3,2,1])));
pz = params.pixel_size.z*1e3;
pxsz = params.pixel_size.x*1e3;
zT = params.fpi.modulation_period*1e3;
ccz = size(coeff,3)/2;
sr.Gainpath = offsetfile;
sr.Peakthresh = 6;
sr.LocTransform = Tm;
sr.Boxsize = bxsz;
sr.Quadshift = trans;
sr.Gain = 2.27;
sr.Phi0 = [0 ,0, 0, 0];
sr.Dz = [0,0,0,0];
sr.Dphi = [0,0,0,0];
sr.BGoffset = 0;
sr.InitPhase = [0,1/3,2/3].*pi;
sr.Initz = [0]*150/pz+ccz;
sr.Initx = [];
sr.Zcenter = ccz;
sr.IABall = coeff;
sr.Pixelsize = pxsz;%nm
pt = zT/2/pi;
sr.Pixelsizez = pz; % nm
sr.Pixelsizeang = pt;% radian to nm
sr.Periodz = zT; % nm
sr.updateTflag = 0;