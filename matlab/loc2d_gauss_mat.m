% for v14 iPALMast_analysisv14scmos.m

function [Res] = loc2d_gauss_mat(obj,datafiles,resfolder)
scmos_cali_file = obj.Gainpath;
gain = obj.Gain;
thresh = obj.Peakthresh;
llthresh = 800;
Ithresh = 100;
%% parameters
warning('off','all')
fn = length(datafiles);
load(scmos_cali_file)
    ccdoffset = mean(single(sequence),3);
    

Res = struct('x',[],'y',[],'t',[],'LL',[],'I',[],'bg',[],'sx',[],'sy',[],'f',[],'stdx',[],'stdy',[],'crlb',[]);
%% load one quadrant of the stack
for ff = 1:fn
    %close all
    disp(['process file number:',num2str(ff)])
    tic
    if iscell(datafiles)
        filename = datafiles{ff};
    else
        filename = datafiles(ff).name;
        datapath = datafiles(ff).folder;
    end
   
filepath = fullfile(datapath,filename);
    F = load(filepath);
    namei = fieldnames(F);
    im1 = (single(F.sequence)-ccdoffset)*gain;
    %im1 = permute(im1,[2,1,3]);
    %% find local centers
    imsz = size(im1);

    if ff == 1
        h = dipshow(max(im1,[],3));
        diptruesize(h,400)
        colormap(grey)
        print(gcf,'-dpng','-r300',[resfolder,obj.Savename,'_sumprojection_image'])
    end

    subsz = obj.Boxsize;   % sub region size


    [rois,cor] = makerois(im1,thresh,subsz);
    numlocs = size(rois,3);
    disp(['A total of ' num2str(numlocs) ' subregions were detected. Start sCMOS_sigmaxy fitting']);

    rois(rois<0) = min(abs(rois(:)));
    

 
    %% find local centers

    
        PSFsigma = 1.6;
    Iterations = 100;
    FitType = 4;
    
    N = size(rois,3);
    maxN = 10000;
    v0 = [1:maxN:N];
    Ns = numel(v0);
    if v0(Ns)<N
        v0(Ns+1) = N;
    end
    v0(end) = v0(end)+1;
    
    P = [];
    CRLB = [];
    LL = [];
    
    for ii = 1:numel(v0)-1
        vi = v0(ii):v0(ii+1)-1;
        [P1, CRLB1, LL1] = GPUgaussMLEv2(single(rois(:,:,vi)),PSFsigma,Iterations,FitType);
        
        P = cat(1,P,P1);
        CRLB = cat(1,CRLB,CRLB1);
        LL = cat(1,LL,LL1);
    end
    
    
    res = struct('x',P(:,2),'y',P(:,1),'I',P(:,3),'bg',P(:,4),'LL',-2*LL,'sx',P(:,5),'sy',P(:,6),...
        'cor',cor,'crlb',CRLB);
    %% filter
    maskxy = res.x<subsz-1 & res.x>1 & res.y<subsz-1 & res.y>1;
    masks = res.sx<subsz/2 & res.sy<subsz/2;
    maskll = res.LL<llthresh; 
    %maskothers = res.I>Ithresh;
    
    mask = maskxy&masks&maskll;
    res1 = applymask(res,mask);
    
    
    %% collect data
    pixelsize = obj.Pixelsize*1e3;
    res1.x = (res1.x+res1.cor(:,2)).*pixelsize;
    res1.y = (res1.y+res1.cor(:,1)).*pixelsize;
    res1.sx = res1.sx.*pixelsize;
    res1.sy = res1.sy.*pixelsize;
    res1.stdx = sqrt(res1.crlb(:,2)).*pixelsize;
    res1.stdy = sqrt(res1.crlb(:,1)).*pixelsize;
    res1.t = res1.cor(:,3);
    res1.f = ones(size(res1.x)).*ff;
    
    Res = catstruct(Res,res1,1);
    toc
end




