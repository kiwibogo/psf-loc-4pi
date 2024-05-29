% for v14 iPALMast_analysisv14scmos.m

function [Res] = loc4pi_psfinv_mac(datapath,datafiles,resfolder,obj)
scmos_cali_file = obj.Gainpath;
thresh = obj.Peakthresh;
gain = obj.Gain;
%% parameters
Res = struct('x',[],'y',[],'zast',[],'phi',[],'t',[],'LL',[],'I',[],'stepN',[],'cycN',[],'bg',[],'sx',[],'sy',[],'ct',[],'cor',[],'convg',[],...
    'f',[],'rms',[],'rmp',[],'stdx',[],'stdy',[],'stdzast',[],'stdphi',[]);
%% load one quadrant of the stack
fn = length(datafiles);

for ff = 1:fn
    %close all
    disp(['process file number:',num2str(ff)])
    tic
    if iscell(datafiles)
        filename = datafiles{ff};
    else
        filename = datafiles(ff).name;
    end
    
    load(fullfile(datapath,filename));
    load(scmos_cali_file);
    
    for ss =1:4
        eval(['ims = double(qd',num2str(ss),')-o',num2str(ss),';']);
        %eval(['ims = double(qd',num2str(ss),')-o(:,:,ss);']);
        ims = ims./gain;
        eval(['qd',num2str(ss),'=ims;']);
    end
    
%     dipshow(cat(1,cat(2,qd1,qd2),cat(2,qd3,qd4)));
    %% find local centers
    imsz = size(qd1);
    shiftxy = obj.Quadshift;

    if ff == 1
        h = figure;
        imagesc(max(qd1,[],3));
        colormap(grey)
        print(gcf,'-dpng','-r300',[resfolder,filename(1:end-4),'_sumprojection_image'])
    end
    
    bxsz = obj.Boxsize;   % sub region size
    T = obj.LocTransform;
    
    [rois,cor,roisum] = makerois_4pi_mac(qd1,qd2,qd3,qd4,T,thresh,bxsz,shiftxy);    
    numlocs = size(rois,3);
    disp(['A total of ' num2str(numlocs) ' subregions were detected. Start sCMOS_sigmaxy fitting']);

    % photometry method 
    PSFsigma = 1.15;
    iterations = 100;
    fittype = 4;
    [P1] = CPUmleFit_LM(single(roisum),fittype,iterations,1,0,0,[PSFsigma,PSFsigma]);
    
    res = struct('x',P1(:,1),'y',P1(:,2),'I',P1(:,3),'bg',P1(:,4),'sx',P1(:,6),'sy',P1(:,5),'cor',cor);
    dx = squeeze(cor(:,1,:)-round(cor(:,1,:)));
    dy = squeeze(cor(:,2,:)-round(cor(:,2,:)));
    
    [rms,rmp,ang_ctr] = getcontrast_mac(qd1,qd2,qd3,qd4,bxsz,res);

    % PSF localiation
    shared = obj.SharedParams;
    rois(rois<0) = min(abs(rois(:)));
    
    [P,CRLB,LogL,PSF] = psfloc_mac(obj,rois,dx,dy,shared);

    
    Ps = cell(1,6);
    CRLBs = cell(1,6);
    n = 0;
    for kk = 1:6
        if shared(kk)==1
            Ps{kk} = P(:,n+kk);
            CRLBs{kk} = CRLB(:,n+kk);
        else
            Ps{kk} = P(:,n+kk:n+kk+3);
            CRLBs{kk} = CRLB(:,n+kk:n+kk+3);
            n = n+3;
        end
    end
    
    res = struct('x',Ps{1},'y',Ps{2},'I',Ps{3},'bg',Ps{4},'LL',-2*LogL,'zast',Ps{5},'phi',Ps{6},'convg',P(:,end),'sx',P1(:,6),'sy',P1(:,5),...
        'cor',cor(:,:,1),'stdx',sqrt(CRLBs{1}),'stdy',sqrt(CRLBs{2}),'stdzast',sqrt(CRLBs{5}),'stdphi',sqrt(CRLBs{6}),'rms',rms,'rmp',rmp,'ct',ang_ctr);
    %% filter
    maskxy = res.x<(bxsz-1)/2+bxsz/4 & res.x>(bxsz-1)/2-bxsz/4 & res.y<(bxsz-1)/2+bxsz/4 & res.y>(bxsz-1)/2-bxsz/4;
    maskxy = mean(maskxy,2)==1;
    maskz = mean(res.zast>1 & res.zast<size(obj.IABall,3)-1,2)==1;
    maskct = res.ct<1.1;
    mask = maskxy & maskct & maskz & res.convg <100;
    res1 = applymask(res,mask);
    res1.stepN = ones(numel(res1.LL),1).*str2double(filename(end-8:end-4))+1;
    res1.cycN = ones(numel(res1.LL),1).*str2double(filename(end-14:end-10))+1;
    
    
    %% collect data
    res1.x = (res1.x+res1.cor(:,1));
    res1.y = (res1.y+res1.cor(:,2));
    res1.zast = (res1.zast-obj.Zcenter);

    res1.t = res1.cor(:,3);
    res1.f = ones(size(res1.LL)).*ff;
    
    Res = catstruct(Res,res1,1);
    
    
    toc
end




