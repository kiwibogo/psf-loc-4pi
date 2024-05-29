% for v14 iPALMast_analysisv14scmos.m

function [Res, obj] = locsingle_psfinv(datapath,datafiles,resfolder,obj)
thresh = obj.Peakthresh;
gain = obj.Gain;
ccdoffset = obj.ccdoffset;
%% parameters

    Res = struct('x',[],'y',[],'z',[],'t',[],'LL',[],'I',[],'bg',[],'cor',[],'convg',[],...
        'f',[],'stdx',[],'stdy',[],'stdz',[]);
%% load one quadrant of the stack
fn = length(datafiles);

for ff = 1:fn
    tic
    %close all

    if iscell(datafiles)
        filename = datafiles{ff};
    else
        filename = datafiles(ff).name;
    end

    filepath = fullfile(datapath,filename);
    F = load(filepath);
    namei = fieldnames(F);
    im1 = (F.(namei{1})-ccdoffset)*gain;
    %im1 = permute(im1,[2,1,3]);
    %% find local centers
    imsz = size(im1);

    if ff == 1
        h = dipshow(max(im1,[],3));
        diptruesize(h,400)
        colormap(grey)
        print(gcf,'-dpng','-r300',[resfolder,obj.Savename,'_sumprojection_image'])
    end

    bxsz = obj.Boxsize;   % sub region size


    [rois,cor] = makerois(im1,thresh,bxsz);
    numlocs = size(rois,3);
    disp(['A total of ' num2str(numlocs) ' subregions were detected. Start sCMOS_sigmaxy fitting']);

    rois(rois<0) = min(abs(rois(:)));
    [P,CRLB,LogL,PSF] = psfloc_ast_single(obj,rois);



    res = struct('x',P(:,2),'y',P(:,1),'I',P(:,3),'bg',P(:,4),'LL',-2*LogL,'z',P(:,5),'convg',P(:,end),...
        'cor',cor,'stdx',sqrt(CRLB(:,2)),'stdy',sqrt(CRLB(:,1)),'stdz',sqrt(CRLB(:,5)));
    %% filter
    maskxy = res.x<(bxsz-1)/2+bxsz/4 & res.x>(bxsz-1)/2-bxsz/4 & res.y<(bxsz-1)/2+bxsz/4 & res.y>(bxsz-1)/2-bxsz/4;
    maskxy = mean(maskxy,2)==1;
    maskz = mean(res.z>1 & res.z<size(obj.coeff,3)-1,2)==1;
    mask = maskxy & maskz & res.convg <max(res.convg);
    res1 = applymask(res,mask);


    %% collect data
    res1.x = (res1.x+res1.cor(:,1));
    res1.y = (res1.y+res1.cor(:,2));
    res1.z = (res1.z-obj.Zcenter);

    res1.t = res1.cor(:,3);
    res1.f = ones(size(res1.LL)).*ff;

    Res = catstruct(Res,res1,1);

    toc
end




