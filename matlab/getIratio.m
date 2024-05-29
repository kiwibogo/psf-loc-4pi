function [sr] = getIratio(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.mat']);
sr.SharedParams = [1,1,0,0,1,1];
%Ns = round(linspace(1,length(datafiles),10));
[res] = loc4pi_psfinv_v1(datapath,datafiles,resultfolder,sr);

maxf = max(res.t)+1;
res.frames = res.t+(res.f-1).*maxf;
sr.Result.init = res;

dt = 0.025;
figure;
subplot(121)
Iavg = movmean(res.I,1e4);
Iratio = Iavg./Iavg(:,1);
plot(res.frames*dt,Iratio)
xlabel('time (s)');
ylabel('I ratio')
legend('1','2','3','4')
subplot(122)
Isp = (Iavg(:,2)+Iavg(:,4))./(Iavg(:,1)+Iavg(:,3));
plot(res.frames*dt,Isp)
xlabel('time (s)');
ylabel('s/p ratio')

sr.Iratio = mean(Iratio);


save([resultfolder, sr.Savename, '_loc_Iratio.mat'],'sr');

%%
figure;
for ss = 2:4
histogram(movmean((res.I(:,ss))./(res.I(:,1)),1),[0:0.05:3]);hold on
end








