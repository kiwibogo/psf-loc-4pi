clearvars
addpath('matlab\')
addpath('mex\')


%% localization

datapath = 'Y:\Projects\SPT-RESOLFT\Data\2022-09-09 - Gamma_mEOS4b\';
resultfolder = datapath;

dataname = 'Cell04_405_0p6V_488_0V_561_10mW_2022_09_09_19_06_02';
fn = length(dataname); 
datafile = [datapath,dataname];
sr.Gainpath = 'Y:\sCMOS Calibrations\SPT\GainCalibration_darkFrames_2022_10_26_18_02_10.mat';
sr.Gain = 0.47;
sr.Boxsize = 11;
 sr.Savename = dataname;
sr.Pixelsize = 0.1; % um
sr.Peakthresh = 1.5;
[res] = fit2d(sr,datafile,resultfolder);

%%
a = 4;
xscale = sr.Pixelsize*1e3/a;
pos = [res.x,res.y,res.t]./[xscale,xscale,1];
dis = pdist(pos);
link = linkage(dis,'single');
Tc = cluster(link,'cutoff',10,'criterion','distance');

clear posM IM
count = 1;
figure;
ha = axes;hold on
for tt=1:max(Tc)
    maskT=(Tc==tt);
    if sum(maskT)>20

        posM{count} = pos(maskT,:)./[a,a,1];
        IM{count} = res.I(maskT);
        plot3(posM{count}(:,1),posM{count}(:,2),posM{count}(:,3),'.-')
        count = count+1;
    end
end




%%
Nt = numel(IM);
err = zeros(1,Nt);
Imin = zeros(1,Nt);
h = figure('Position',[200,100,600,500]);
tg1 = uitabgroup(h); % tabgroup
for ii=1:Nt
    trk = round(posM{ii});
    I1 = double(IM{ii});
    I1(I1>4*median(I1)) = min(I1);
    I = I1-min(I1);
    I = I./max(I);
    Iavg = movmean(I,3);
    startpoint = [0.7,0.3,0.5];
    Ns = 2;
    est = fminsearch(@(x)stepfit(x,Ns,I),startpoint);

    [mse,Ifit] = stepfit(est,Ns,I);
    err(ii) = mse;
    Imin(ii) = est(1)*(max(I1)-min(I1))+min(I1);
    tab1 = uitab(tg1,'title',num2str(ii));
    ha = axes('parent',tab1);
    plot(trk(:,3),[I,Ifit'])
end
%%
ind = [1,19,21,32,33,37,41];
figure;
plot(err,'.-');hold on;
plot(ind,err(ind),'o')

%%
icp = smi_stat.ChangeDetection(round(double(IM{1})),200);
icp.plotIntensityEstimate()
%%
load(sr.Gainpath)
ccdoffset = mean(single(sequence),3);
F = load([datafile,'.mat']);
imsz = size(sequence);
im = (single(F.sequence)-ccdoffset);
%%
ind = 1;
trk = round(posM{ind});
imt = zeros(imsz);
for ii = 1:size(trk,1)
imt(trk(ii,1),trk(ii,2),trk(ii,3)+1) = 1;
end
overlay(im.*2,imt)



%%
function [mse,I] = stepfit(x,Ns,data)
I0 = x(1:Ns);
ts = x(Ns+1:2*Ns-1);
L = length(data);
I = fstep(I0,Ns,L,ts);
mse = mean((I(:)-data(:)).^2)/mean(data);
end

function I = fstep(I0,Ns,L,ts)

t = [0:L-1];
I = zeros(1,L);

ts(ts<0.2) = 0.2;
ts(ts>0.8) = 0.8;
ind = [1,round(ts.*L),L];

for ii = 1:Ns
    I(ind(ii):ind(ii+1)) = (Ns-ii+1)*I0(1);
end

end


