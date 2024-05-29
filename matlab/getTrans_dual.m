function [sr] = getTrans_dual(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.tif']);
sr.Iratio = [1,1];
sr.SharedParams = [0,0,1,0,0];
sr.updateT = 1;

[res, sr] = locdual_psfinv(datapath,datafiles(1),resultfolder,sr,[2]);

dx = res.x-res.x(:,1);
dy = res.y-res.y(:,1);
figure;
mask = dx(:,2)>quantile(dx(:,2),0.1) & dx(:,2)<quantile(dx(:,2),0.9) & dy(:,2)>quantile(dy(:,2),0.1) & dy(:,2)<quantile(dy(:,2),0.9);
subplot(121);plot(movmean(dx(mask,:),1000));xlabel('localizations');ylabel('dx')
subplot(122);plot(movmean(dy(mask,:),1000));xlabel('localizations');ylabel('dy')
legend('1','2')

sr.DataShift = [median(dx,1);median(dy,1)];

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
