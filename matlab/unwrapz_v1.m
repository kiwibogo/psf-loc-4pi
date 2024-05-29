function [sr,res] = unwrapz_v1(sr,dt,resultfolder)
res = sr.Result.init;
res.zast = res.zast*sr.Pixelsizez;
zastlim = [-400,400];
mask = res.ct>0.3 & res.zast>zastlim(1) & res.zast<zastlim(2);
res1 = applymask(res,mask);
zT = sr.Periodz;
w = 6000;
df = 4000;
N = length(res1.zast);
Nf = length(1:df:N-w);
z0all = zeros(Nf+2,1);
fcc = zeros(Nf+2,1);
phioffset = 0;
count = 1;
for ii = [0,1:df:N-w,N]
    if ii==0
        ind = [1:round(w/1.5)];
        fcc(count) = 0;
    elseif ii==N
        ind = [N-round(w/1.5):N];
        fcc(count) = res1.frames(N);
    else
        ind = [ii:ii+w];
        fcc(count) = mean(res1.frames(ind));
    end
    
    zast = res1.zast(ind);
    phi = mod(res1.phi(ind),2*pi);
    
    zfp = (phi+phioffset)/2/pi*zT;
    dz = zfp-zast;
    dzm = mod(dz,zT);
    z0 = -cyclicaverage(dzm,zT)+phioffset/2/pi*zT;
    z0all(count) = z0;
    count = count+1;
end

z0uw = unwrap(z0all./zT.*2*pi)./2/pi*zT;

fz = fit(fcc,z0uw,'smoothingspline','SmoothingParam',0.5e-8);
z0fit = fz(res.frames);
figure;plot(fcc*dt,z0uw,'.-',res.frames*dt,z0fit); ylabel('phase drift (nm)'); xlabel('time (s)')
%%
phi = mod(res.phi,2*pi);
phical = (res.zast-z0fit)./zT*2*pi;
dphi = -(phi-phical)./pi/2;
period = round((dphi));
res.zphi = period.*zT+phi./2/pi*zT;


fnum = 50;
[cormask,fcc] = gencormask(res.frames,res.frames,fnum);

z0avg = zeros(size(z0fit));
for ii = 1:max(cormask)
    mask = cormask == ii;
    z0avg(mask) = mean(z0fit(mask));    
end

figure;plot(res.frames,z0avg,'.-')

res.z = res.zphi+z0avg;

%%
res.x = res.x.*sr.Pixelsize;
res.y = res.y.*sr.Pixelsize;
res.stdx = res.stdx.*sr.Pixelsize;
res.stdy = res.stdy.*sr.Pixelsize;
res.stdzast = res.stdzast*sr.Pixelsizez;

res.phi0 = z0fit./zT.*2*pi;
res.stdz = res.stdphi.*sr.Pixelsizeang;
res.stdxy = sqrt(res.stdx.^2+res.stdy.^2);

%% get density map
fnum = 10;
[cormask,fcc] = gencormask(res.zast-min(res.zast),res.zast-min(res.zast),fnum);
fcc = fcc+min(res.zast);
Ns = length(fcc);
figure('position',[200,200,600,500]);
ha = dscatter(res.zast,res.z);
density = ha.Children.CData;
zpeak = zeros(Ns,1);
id = unique(cormask);
for ii = 1:Ns
    mask = cormask==id(ii);
    z = res.z(mask);
    ds = density(mask);
    zpeak(ii) = mean(z(find(ds==max(ds))));
end

hold on; plot(fcc,zpeak,'r.-');xlabel('zastig (nm)');ylabel('z phase (nm)')

fzp = fit(fcc,zpeak,'smoothingspline','SmoothingParam',1e-4);
zpfit = fzp(res.zast);
res.zerr = abs(res.z-zpfit);
res.density = density;
sr.Result.unwrap = res;

save([resultfolder, sr.Savename, '_loc_z.mat'],'res','sr');
save([resultfolder, sr.Savename, '_loc_z_smap.mat'],'res');


