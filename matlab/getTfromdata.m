function [sr] = getTfromdata(sr,datapath,resultfolder,fitflag)
datafiles = dir([datapath,'*.mat']);
if fitflag == 1
    sr.SharedParams = [0,0,1,1,1,1];
    [res] = loc4pi_psfinv_v1(datapath,datafiles,resultfolder,sr);
    maxf = max(res.t)+1;
    res.frames = res.t+(res.f-1).*maxf;
    sr.Result.trans = res;    
    save([resultfolder, sr.Savename, '_loc_ini_trans.mat'],'sr');
end
dx = movmean(res.x-res.x(:,1),2400);
dy = movmean(res.y-res.y(:,1),2400);
dt = 0.025; % s
pxsz = 129;
figure;
subplot(121);plot(res.frames*dt,dx.*pxsz,'-');xlabel('localizations');ylabel('dx')
subplot(122);plot(res.frames*dt,dy.*pxsz,'-');xlabel('localizations');ylabel('dy')
legend('1','2','3','4')

N1 = length(res.frames);
dxfit = zeros(N1,4);
dyfit = zeros(N1,4);
fitrange = round(linspace(1,N1,1000));
T = struct('x',[],'y',[]);
for ss = 1:4
    fx = fit(res.frames(fitrange),dx(fitrange,ss),'smoothingspline','SmoothingParam',1e-9);
    T.x{ss} = fx;
    dxfit(:,ss) = fx(res.frames);
    fy = fit(res.frames(fitrange),dy(fitrange,ss),'smoothingspline','SmoothingParam',1e-9);
    T.y{ss} = fy;
    dyfit(:,ss) = fy(res.frames);
end


figure('position',[200,300,900,400]);
ha = subplot(121);plot(res.frames*dt,dx.*pxsz,'-');hold on;plot(res.frames*dt,dxfit.*pxsz,'k--');ylabel('dx (nm)');xlabel('time (s)')
ha = subplot(122);plot(res.frames*dt,dy.*pxsz,'-');hold on;plot(res.frames*dt,dyfit.*pxsz,'k--');ylabel('dy (nm)');xlabel('time (s)')

legend('1','2','3','4')

sr.Transoffset = T;

save([resultfolder, sr.Savename, '_loc_ini_trans.mat'],'sr');
