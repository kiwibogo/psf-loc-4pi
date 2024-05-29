function dcimg2mat(datapath,datafiles,centers,singlefile)
fn = size(datafiles,1);
imsL = [];
qdL1 = [];
qdL2 = [];
qdL3 = [];
qdL4 = [];
try
    metadata = ini2struct([datapath,'\','metadata.ini']);
    metadata.microscopeID = '4Pi';
catch
    metadata.microscopeID = '4Pi';
end

if singlefile == 0
    for ii =1:fn
        filename = datafiles(ii).name;
        if isempty(centers)
            ims = iPALM_readdcimg1(fullfile(datapath,filename));
            save([fullfile(datapath,filename(1:end-6)),'.mat'],'ims')
        else
            [qds]=iPALM_readdcimg1(fullfile(datapath,filename),centers);
            qd1 = qds(:,:,:,1);
            if numel(centers(1,:))==4
                qd2 = qds(:,:,:,2);
                qd3 = qds(:,:,:,3);
                qd4 = qds(:,:,:,4);
                %filename(end-13:end-12) = num2str(ii+89);
                save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','qd2','qd3','qd4','metadata')
            elseif numel(centers)==2
                save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','metadata')
            else
                qd2 = qds(:,:,:,2);
                save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','qd2','metadata')
            end
        end
        clear qd1 qd2 qd3 qd4 qds
        delete(fullfile(datapath,filename))
    end
    
else
    if isempty(centers)
        for ii = 1:fn
            filename = datafiles(ii).name;
            ims1 = iPALM_readdcimg1(fullfile(datapath,filename));
            if ii>1
                sz1 = size(ims1);
                sz2 = size(imsL);
                if numel(sz1)<3 || numel(sz2)<3
                    sz2(3) = 1;
                    sz1(3) = 1;
                end
                ims1 = ims1(:,:,1:min([sz1(3),sz2(3)]),:);
                sz3 = size(ims1);
                if numel(sz3)<3
                    sz3(3)=1;
                end
                imsL = cat(4,imsL(:,:,1:sz3(3),:),ims1);
            else
                imsL = cat(4,imsL,ims1);
            end
            delete(fullfile(datapath,filename))
        end
        ims = imsL;
        save([fullfile(datapath,filename(1:end-6)),'.mat'],'ims','metadata')
    else
        for ii = 1:fn
            filename = datafiles(ii).name;
            [qds]=iPALM_readdcimg1(fullfile(datapath,filename),centers);
            if ii>1
                sz1 = size(qds);
                sz2 = size(qdL1);
                if numel(sz1)<3 || numel(sz2)<3
                    sz1(3) = 1;
                    sz2(3) = 1;
                end
                qds = qds(:,:,1:min([sz1(3),sz2(3)]),:);
                qdL1 = qdL1(:,:,1:min([sz1(3),sz2(3)]),:);
                qdL2 = qdL2(:,:,1:min([sz1(3),sz2(3)]),:);
                qdL3 = qdL3(:,:,1:min([sz1(3),sz2(3)]),:);
                qdL4 = qdL4(:,:,1:min([sz1(3),sz2(3)]),:);
            end
            qd1 = qds(:,:,:,1);
            qdL1 = cat(4,qdL1,qd1);
            if numel(centers(1,:))==4
                qd2 = qds(:,:,:,2);
                qd3 = qds(:,:,:,3);
                qd4 = qds(:,:,:,4);
                qdL2 = cat(4,qdL2,qd2);
                qdL3 = cat(4,qdL3,qd3);
                qdL4 = cat(4,qdL4,qd4);
                elseif numel(centers(1,:))==2 && numel(centers)==4
                    qd2 = qds(:,:,:,2);
                    qdL2 = cat(4,qdL2,qd2);
            end
            delete(fullfile(datapath,filename))
        end
        qd1 = qdL1;
        
        if numel(centers(1,:))==4
            qd2 = qdL2;
            qd3 = qdL3;
            qd4 = qdL4;
            save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','qd2','qd3','qd4','metadata')
        elseif numel(centers)==2
            save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','metadata')
        else
            save([fullfile(datapath,filename(1:end-6)),'.mat'],'qd1','qd2','metadata')
            
        end
    
        
    end
end