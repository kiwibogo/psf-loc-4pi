function res = fit2d(obj,dataname,Resultpath)
resfolder = [Resultpath,'\2D fit\'];
if ~exist(resfolder,'dir')
mkdir(resfolder)
end
%%
datafiles = dir([dataname,'*.mat']);
[res] = loc2d_gauss_mat(obj,datafiles,resfolder);
maxf = max(res.t)+1;
res.frames = res.t+(res.f-1).*maxf;
obj.Savename = datafiles(1).name(1:end-6);
obj.Result2d.init = res;
save([resfolder, datafiles(1).name(1:end-6), '_loc_2d'],'res');

end