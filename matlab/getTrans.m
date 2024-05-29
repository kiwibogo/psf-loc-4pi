function [sr] = getTrans(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.mat']);
sr.Peakthresh = 20;
sr.Iratio = [1,1,1,1];
sr.SharedParams = [0,0,1,1,1,1];
[res] = loc4pi_psfinv_v1(datapath,datafiles(2),resultfolder,sr);

dx = res.x-res.x(:,1);
dy = res.y-res.y(:,1);
mask = abs(dx(:,1))>-1;
for ss = 2:4
mask = mask & dx(:,ss)>quantile(dx(:,ss),0.1) & dx(:,ss)<quantile(dx(:,ss),0.9) & dy(:,ss)>quantile(dy(:,ss),0.1) & dy(:,ss)<quantile(dy(:,ss),0.9);
end
figure;
subplot(121);plot(movmean(dx(mask,:),1000));xlabel('localizations');ylabel('dx')
subplot(122);plot(movmean(dy(mask,:),1000));xlabel('localizations');ylabel('dy')
legend('1','2')

sr.DataShift = [median(dy,1);median(dx,1)];


%%
% maxf = max(res.t)+1;
% res.frames = res.t+(res.f-1).*maxf;
% 
% dt = 0.025;
% t = res.frames*dt;
% h = figure('position',[200,300,450,400]);
% plot(t,Iavg./Iavg(:,1),'linewidth',1.2);
% xlim([0,30])
% xlabel('time (s)');ylabel('intensity ratio')
% legend('1','2','3','4')
% 
% set(h,'PaperPositionMode','auto')
% imgname = 'Iratio';
% print(h,'-depsc','-r300',[resultfolder,imgname])
% print(h,'-dpng','-r300',[resultfolder,imgname])
