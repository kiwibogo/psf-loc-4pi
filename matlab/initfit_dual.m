function [sr,res] = initfit_dual(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.tif']);

sr.SharedParams = [1,1,1,0,0];
sr.updateT = 0;
[res] = locdual_psfinv(datapath,datafiles,resultfolder,sr);

%%
maxf = max(res.t)+1;
res.frames = res.t+(res.f-1).*maxf;
res.x = res.x.*sr.Pixelsizex;
res.y = res.y.*sr.Pixelsizey;
res.z = res.z.*sr.Pixelsizez;

res.stdx = res.stdx.*sr.Pixelsizex;
res.stdy = res.stdy.*sr.Pixelsizey;
res.stdz = res.stdz.*sr.Pixelsizez;


sr.Result.init = res;

save([resultfolder, sr.Savename, '_loc.mat'],'sr');

