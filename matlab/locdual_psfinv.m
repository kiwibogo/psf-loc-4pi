% for v14 iPALMast_analysisv14scmos.m

function [Res, obj] = locdual_psfinv(datapath,datafiles,resfolder,obj,fitrange)
thresh = obj.Peakthresh;
gain = obj.Gain;
ccdoffset = obj.ccdoffset;
T = obj.LocTransform;
%% parameters
Res = struct('x',[],'y',[],'z',[],'t',[],'LL',[],'I',[],'bg',[],'cor',[],'convg',[],...
    'f',[],'stdx',[],'stdy',[],'stdz',[]);
%% load one quadrant of the stack
fn = length(datafiles);

for ff = 1:fn
    %close all
    
    if iscell(datafiles)
        filename = datafiles{ff};
    else
        filename = datafiles(ff).name;
    end
    
    filepath = fullfile(datapath,filename);

    tiff_info = imfinfo(filepath); % return tiff structure, one element per image
    xsz = tiff_info(1).Width;
    ysz = tiff_info(1).Height;
    fmax = length(tiff_info);
    Nf = 2000;
    fID = fopen (filepath, 'r');
    if nargin<5
        fitrange = 1:ceil(fmax/Nf);
        if ff == 1
            fitrange([1]) = [];
        end
    end

    for nn = fitrange
        tic
        fstart = (nn-1)*Nf+1;
        disp(['process file number ',num2str(ff),' frame ', num2str(fstart)])
        ims = zeros(ysz,xsz,Nf,'uint16');
        for ii = fstart : min([fstart+Nf,fmax])
            fseek (fID, tiff_info(ii).StripOffsets, 'bof');
            %ims(:,:,ii) = imread(filepath, ii);
            ims(:,:,ii-fstart+1) = fread(fID,[xsz,ysz],'uint16=>uint16')';
        end
        ims = (single(ims)-ccdoffset)/gain;
        midpoint = round(ysz/2);
        im2 = flip(ims(1:midpoint,:,:),1);
        im1 = ims(midpoint+1:end,:,:); %reference channel is im1
               
        im1 = flip(im1,2);
        im2 = flip(im2,2);
        %% find local centers
        imsz = size(im1);
        
        if ff == 1 && nn==2
            h = dipshow(max(im1,[],3));
            diptruesize(h,400)
            colormap(grey)
            print(gcf,'-dpng','-r300',[resfolder,obj.Savename,'_sumprojection_image'])
        end
        
        bxsz = obj.Boxsize;   % sub region size
        

        [rois,cor,T] = makerois_dual(im1,im2,T,thresh,bxsz, obj.updateT);
        numlocs = size(rois,3);
        disp(['A total of ' num2str(numlocs) ' subregions were detected. Start sCMOS_sigmaxy fitting']);

        
        dx = squeeze(cor(:,1,:)-round(cor(:,1,:)));
        dy = squeeze(cor(:,2,:)-round(cor(:,2,:)));
        

        
        
        shared = obj.SharedParams;
        rois(rois<0) = min(abs(rois(:)));

        [P,CRLB,LogL] = psfloc_dual(obj,rois,dy,dx,shared);
        
        
        Ps = cell(1,5);
        CRLBs = cell(1,5);
        n = 0;
        for kk = 1:5
            if shared(kk)==1
                Ps{kk} = P(:,n+kk);
                CRLBs{kk} = CRLB(:,n+kk);
            else
                Ps{kk} = P(:,n+kk:n+kk+1);
                CRLBs{kk} = CRLB(:,n+kk:n+kk+1);
                n = n+1;
            end
        end
        
        res = struct('x',Ps{2},'y',Ps{1},'I',Ps{4},'bg',Ps{5},'LL',-2*LogL,'z',Ps{3},'convg',P(:,end),...
            'cor',cor(:,:,1),'stdx',sqrt(CRLBs{1}),'stdy',sqrt(CRLBs{2}),'stdz',sqrt(CRLBs{3}));
        %% filter
        maskxy = res.x<(bxsz-1)/2+bxsz/4 & res.x>(bxsz-1)/2-bxsz/4 & res.y<(bxsz-1)/2+bxsz/4 & res.y>(bxsz-1)/2-bxsz/4;
        maskxy = mean(maskxy,2)==1;
        maskz = mean(res.z>1 & res.z<size(obj.coeff,3)-1,2)==1;
        mask = maskxy & maskz & res.convg <max(res.convg)-5;
        res1 = applymask(res,mask);
        
        
        %% collect data
        res1.x = (res1.x+res1.cor(:,1));
        res1.y = (res1.y+res1.cor(:,2));
        res1.z = (res1.z-obj.Zcenter);
        
        res1.t = res1.cor(:,3)+fstart;
        res1.f = ones(size(res1.LL)).*ff;
        
        Res = catstruct(Res,res1,1);
        
        
        toc
    end
end
obj.LocTransform = T;




