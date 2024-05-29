function [sr,res] = unwrapz_v2(sr,res,dt)
zastlim = [-450,450];
mask = res.zastig>zastlim(1) & res.zastig<zastlim(2) & res.iterations<100;
res1 = applymask(res,mask);
zT = sr.Periodz;
w = 6000;
df = 4000;
N = length(res1.zastig);
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
        fcc(count) = res1.frame(N);
    else
        ind = [ii:ii+w];
        fcc(count) = mean(res1.frame(ind));
    end
    
    zast = res1.zastig(ind);
    phi = mod(res1.phase(ind),2*pi);
    
    zfp = (phi+phioffset)/2/pi*zT;
    dz = zfp-zast;
    dzm = mod(dz,zT);
    z0 = -cyclicaverage(dzm,zT)+phioffset/2/pi*zT;
    z0all(count) = z0;
    count = count+1;
end

z0uw = unwrap(z0all./zT.*2*pi)./2/pi*zT;

fz = fit(fcc,z0uw,'smoothingspline','SmoothingParam',0.5e-8);
z0fit = fz(res.frame);
figure;plot(fcc*dt,z0uw,'.-',res.frame*dt,z0fit); ylabel('phase drift (nm)'); xlabel('time (s)')
%%
phi = mod(res.phase,2*pi);
phical = (res.zastig-z0fit)./zT*2*pi;
dphi = -(phi-phical)./pi/2;
period = round((dphi));
res.zphase = period.*zT+phi./2/pi*zT;


fnum = 50;
[cormask,fcc] = gencormask(res.frame,res.frame,fnum);

z0avg = zeros(size(z0fit));
for ii = 1:max(cormask)
    mask = cormask == ii;
    z0avg(mask) = mean(z0fit(mask));    
end

figure;plot(res.frame,z0avg,'.-')

res.znm = res.zphase+z0avg;

%%

res.phi0 = z0fit./zT.*2*pi;
res.znmerr = res.phaseerr.*sr.Pixelsizeang;

%% get density map
fnum = 10;
[cormask,fcc] = gencormask(res.zastig-min(res.zastig),res.zastig-min(res.zastig),fnum);
fcc = fcc+min(res.zastig);
Ns = length(fcc);
figure('position',[200,200,600,500]);
ha = dscatter(res.zastig,res.znm);
density = ha.Children.CData;
zpeak = zeros(Ns,1);
id = unique(cormask);
for ii = 1:Ns
    mask = cormask==id(ii);
    z = res.znm(mask);
    ds = density(mask);
    zpeak(ii) = mean(z(find(ds==max(ds))));
end

hold on; plot(fcc,zpeak,'r.-');xlabel('zastig (nm)');ylabel('z phase (nm)')

fzp = fit(fcc,zpeak,'smoothingspline','SmoothingParam',1e-4);
zpfit = fzp(res.zastig);
res.zerr = abs(res.znm-zpfit);
res.density = density;
sr.Result.unwrap = res;



