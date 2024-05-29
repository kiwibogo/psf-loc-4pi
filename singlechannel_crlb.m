clearvars
addpath('matlab\')
addpath('mex\')
%addpath('C:\Users\Ries Lab\git\SMAP\shared\myfunctions')

%% load model


mfilepath = 'Y:\Personal Folders\Sheng\data\PSF Engineering Microscope\PR_10-04-2022_UAF_7px_Shift_200nm_beads_Gain0_PupilSeq\';
modelname = 'ZStack_-2to+2um_PupilSize125_psfmodel_pupil_vector_single.h5';
Fm = loadh5([mfilepath,modelname]);
val = h5readatt([mfilepath,modelname],'/','params');
params = jsondecode(val);
pupil = Fm.res.pupil;
%% generate index mismatch model
stagepos = 2; % um
zrange = [-stagepos:0.05:0.5]; % um

%[PSFs,zemit] = vectorPSF_single(Fm,params, zrange, stagepos);
%%
stagepos = 0.5; % um
zrange = [-stagepos:0.05:0.5]; % um
fillfactor = 0.82;
shiftx = -2.5;
pupil_mag = (abs(Fm.res.pupil)>0);
Fm.res.pupil = pupil_mag;
params.option.imaging.NA = 1.7;
params.option.imaging.RI.imm = 1.78;
[PSFs,zemit] = vectorPSF_single(Fm,params, zrange, stagepos, fillfactor,shiftx,'ua');

% generate spline coeff
coeff = Spline3D_interp(PSFs);
%
photon = 50000;
bg = 10;
data = poissrnd(PSFs*photon+bg);
bxsz = size(data,1);
sr = resetparam_imm(params,bxsz,coeff);
[P,CRLB,LogL,PSF] = psfloc_ast_single(sr,data);

stdM = sqrt(CRLB(:,[1,2,5])).*[sr.Pixelsizex,sr.Pixelsizey,sr.Pixelsizez];

figure('position',[200,200,500,400]);
plot(zemit,stdM);
title(['stage position: ',num2str(stagepos),' \mum'])
xlabel('emitter z (um)'); ylabel('std (nm)')
ylim([0,12])
stdM(1,:)

%% vary shift
figure('position',[200,200,500,400]);
ha = axes;hold on;
stagepos = 1; % um
zrange = [-stagepos:0.05:0.5]; % um
fillfactor = 0.846;
shiftx = [-2.0,-2.5,-3,-3.5];
pupil_mag = (abs(Fm.res.pupil)>0);
Fm.res.pupil = pupil_mag;
for ii = 1:numel(shiftx)
[PSFs,zemit] = vectorPSF_single(Fm,params, zrange, stagepos, fillfactor,shiftx(ii),'ua');

% generate spline coeff
coeff = Spline3D_interp(PSFs);
%
photon = 50000;
bg = 10;
data = poissrnd(PSFs*photon+bg);
bxsz = size(data,1);
sr = resetparam_imm(params,bxsz,coeff);
[P,CRLB,LogL,PSF] = psfloc_ast_single(sr,data);

stdM = sqrt(CRLB(:,[1,2,5])).*[sr.Pixelsizex,sr.Pixelsizey,sr.Pixelsizez];


plot(zemit,stdM(:,3));
title(['stage position: ',num2str(stagepos),' \mum'])
xlabel('emitter z (um)'); ylabel('std (nm)')
ylim([1,2])
xlim([0,0.5])
labels{ii} = ['shift ',num2str(shiftx(ii)),' pixels'];
end
legend(labels)

%% vary fillfactor
figure('position',[200,200,500,400]);
ha = axes;hold on;
stagepos = 1; % um
zrange = [-stagepos:0.05:0.5]; % um
fillfactor = [0.83,0.845,0.86,0.865];
shiftx = -2.5;
pupil_mag = (abs(Fm.res.pupil)>0);
Fm.res.pupil = pupil_mag;
for ii = 1:numel(fillfactor)
[PSFs,zemit] = vectorPSF_single(Fm,params, zrange, stagepos, fillfactor(ii),shiftx,'ua');

% generate spline coeff
coeff = Spline3D_interp(PSFs);
%
photon = 50000;
bg = 10;
data = poissrnd(PSFs*photon+bg);
bxsz = size(data,1);
sr = resetparam_imm(params,bxsz,coeff);
[P,CRLB,LogL,PSF] = psfloc_ast_single(sr,data);

stdM = sqrt(CRLB(:,[1,2,5])).*[sr.Pixelsizex,sr.Pixelsizey,sr.Pixelsizez];


plot(zemit,stdM(:,3));
title(['stage position: ',num2str(stagepos),' \mum'])
xlabel('emitter z (um)'); ylabel('std (nm)')
ylim([1,2])
xlim([0,0.5])
labels{ii} = ['fill factor ',num2str(fillfactor(ii))];
end
legend(labels)