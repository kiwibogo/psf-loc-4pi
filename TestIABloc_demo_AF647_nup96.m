clearvars
addpath('matlab\')
addpath('matlab\easyh5\')
addpath('mex\')
addpath('C:\Users\Ries Lab\git\SMAP\shared\myfunctions')
offsetfile = 'test data\ccdoffset.mat';

%% load model
modelpath = 'D:\Sheng\data\06-15-2022 seipin AF647\100nm_676\';
modelname = 'psfmodel_voxel_4pi.h5';
filename = fullfile(modelpath,modelname);
Fm = loadh5(filename);
val = h5readatt(filename,'/','params');
params = jsondecode(val{1});

%% localization

datapath = 'D:\Sheng\data\06-15-2022 seipin AF647\';
resultfolder = [datapath,'\psf fit\'];
if ~exist(resultfolder,'dir')
    mkdir(resultfolder)
end

dataname = {'cell5_seipin_013'};
fn = length(dataname); 

file1 = dir([datapath,dataname{1},'*']);
datapath1 = [datapath,file1(1).name,'/'];
datafiles = dir([datapath1,'*.mat']);
load(fullfile(datapath1,datafiles(1).name));
shiftxy = findshift(max(qd1,[],3),max(qd3,[],3),'iter');
 
bxsz = 13;% subregion size
sr = resetparam_h5(offsetfile,Fm,params,bxsz,shiftxy');
for ii = 1:fn

    file1 = dir([datapath,dataname{ii},'*']);
    datapath1 = [datapath,file1(1).name,'/'];
    datafiles = dir([datapath1,'*.mat']);
    [~,savename] = fileparts(datafiles(1).name);
    sr.Savename = savename;
    % update traslation in Transformation matrix
    [sr] = getTrans(sr,datapath1,resultfolder);
    sr.LocTransform(1:2,3,:) = sr.LocTransform(1:2,3,:)+reshape(sr.DataShift,2,1,4);
        
    
    % initial fit
    sr.Peakthresh = 10;
    [sr] = initfit(sr,datapath1,resultfolder);
    
    % unwrapping
    dt = str2double(metadata.Acquisition.CaptureExposureTime); % exposure time
    [sr] = unwrapz_v1(sr,dt,resultfolder);
    
end
%% dme drift correction
% first run python code to get shift in xyz, saved as *dmeshift.mat

for nn = 1:fn
dmefile = dir([resultfolder,dataname{nn},'*dmeshift.mat']);
load([resultfolder,dmefile(1).name])
resfile = dir([resultfolder,dataname{nn},'*smap.mat']);
load([resultfolder,resfile(1).name])
px = str2double(metadata.PhaseRetrieve.pixelSize)*1e3;
g = [px,px,1000]; % pixel size in x,y,z, unit: nm
shift_dme = shift_dme-shift_dme(1,:);
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
dt = str2double(metadata.Acquisition.CaptureExposureTime); % s
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

