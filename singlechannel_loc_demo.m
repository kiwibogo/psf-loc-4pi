clearvars
addpath('matlab\')
addpath('mex\')
addpath('C:\Users\Ries Lab\git\SMAP\shared\myfunctions')

%% load model


mfilepath = 'Y:\Personal Folders\Sheng\data\PSF Engineering Microscope\PR_10-04-2022_UAF_7px_Shift_200nm_beads_Gain0_PupilSeq\';
modelname = 'ZStack_-2to+2um_PupilSize125_psfmodel_pupil_vector_single.h5';
Fm = loadh5([mfilepath,modelname]);
val = h5readatt([mfilepath,modelname],'/','params');
params = jsondecode(val);

%% generate index mismatch model
stagepos = 1; % um
zrange = [-1:0.05:1]; % um

[PSFs] = vectorPSF_single(Fm,params, zrange, stagepos);


%% generate spline coeff
coeff = Spline3D_interp(PSFs);

%% localization
%datapath = 'E:\Yiming\190910_beads3Ddualcolor_M2\';
datapath = 'Y:\Personal Folders\Sheng\data\PSF Engineering Microscope\PR_10-04-2022_UAF_7px_Shift_200nm_beads_Gain0_PupilSeq\';
resultfolder = datapath;

dataname = {'pupil'};
fn = length(dataname); 

coeff = permute(Fm.locres.coeff,[4,3,2,1]);

bxsz = 20;
sr = resetparam_imm(params,bxsz,coeff);
%sr.coeff = coeff;

for ii = 1:fn
    
    file1 = dir([datapath,dataname{ii},'*']);
    datapath1 = [datapath,file1(1).name,'/'];
    datafiles = dir([datapath1,'*.mat']);
    [~,savename] = fileparts(datafiles(1).name);
    [~,savename] = fileparts(savename);
    sr.Savename = [savename];
    
    
    % initial fit
    sr.Peakthresh = 50;

    
    [sr] = fit_single(sr,datapath1,resultfolder);
    
end




%% dme drift correction
% first run python code to get shift in xyz, saved as *dmeshift.mat
for nn = 1:fn
dmefile = dir([resultfolder,dataname{nn},'*I_smap*dmeshift.mat']);
load([resultfolder,dmefile(1).name])
resfile = dir([resultfolder,dataname{nn},'*I_smap.mat']);
load([resultfolder,resfile(1).name])
g = [127,127,1000]; % pixel size in x,y,z, unit: nm
shiftavg = movmean(shift_dme,200);
fcc = [min(res.frames):max(res.frames)]';


fitrange = round(linspace(1,numel(fcc),1000));
dxyz = [];
for ii = 1:3
fz = fit(fcc(fitrange),shiftavg(fitrange,ii),'smoothingspline','SmoothingParam',1e-7);
dxyz = cat(2,dxyz,fz(fcc));
end

%
label = {'x','y','z'};
dt = 0.025; % s
t1 = fcc*dt./60;
h = figure;
h.Position = [230,300,1200,920];
for ii = 1:3
ha = subplot(3,1,ii);hold on
plot(t1,shift_dme(:,ii).*g(ii),'color',[1,1,1]*0.7,'linewidth',1.2);
plot(t1,dxyz(:,ii).*g(ii),'linewidth',1.2);
ha.Title.String = label{ii};
ha.XLabel.String = 'time (min)';
ha.YLabel.String = 'drift (nm)';
end

res1 = res;
frames = res.frames-min(res.frames)+1;
res1.x = res.x-dxyz(frames,1)*g(1);
res1.y = res.y-dxyz(frames,2)*g(2);
res1.z = res.z-dxyz(frames,3)*g(3);
save([resultfolder,resfile(1).name(1:end-4),'_dme.mat'],'res1')

end

