function [sr,res] = globfit_dual(sr,datapath,resultfolder, photon_ratio)
datafiles = dir([datapath,'*.tif']);

sr.SharedParams = [1,1,1,1,0];
sr.updateT = 0;
sr.photon_ratio = photon_ratio;
[res] = locdual_psfinv_glob(datapath,datafiles,resultfolder,sr);
    
    
    

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


sr.Result.glob = res;

save([resultfolder, sr.Savename, '_loc.mat'],'sr');
save([resultfolder, sr.Savename, '_loc_I_smap.mat'],'res');
