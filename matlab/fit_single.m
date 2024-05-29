function [sr,res] = fit_single(sr,datapath,resultfolder)
datafiles = dir([datapath,'*.mat']);


[res] = locsingle_psfinv(datapath,datafiles,resultfolder,sr);
    
    
    

%%
maxf = max(res.t)+1;
res.frames = res.t+(res.f-1).*maxf;
res.frames = res.frames-min(res.frames);
res.x = res.x.*sr.Pixelsizex;
res.y = res.y.*sr.Pixelsizey;
res.z = res.z.*sr.Pixelsizez;

res.stdx = res.stdx.*sr.Pixelsizex;
res.stdy = res.stdy.*sr.Pixelsizey;
res.stdz = res.stdz.*sr.Pixelsizez;


sr.Result = res;

save([resultfolder, sr.Savename, '_loc.mat'],'sr');
save([resultfolder, sr.Savename, '_loc_I_smap.mat'],'res');
