function [sr,res] = initfit(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.mat']);

sr.SharedParams = [1,1,1,1,1,1];
sr.updateTflag = 0;
[res] = loc4pi_psfinv_v1(datapath,datafiles,resultfolder,sr);

%%
maxf = max(res.t)+1;
res.frames = res.t+(res.f-1).*maxf;
sr.Result.init = res;

save([resultfolder, sr.Savename, '_loc_inifit.mat'],'sr');
