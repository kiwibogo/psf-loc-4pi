function [qds]=iPALMscmos_makeqds(ims,center,fsign)

vsz=size(ims,1);
N = size(ims,3);
qds=uint16(zeros(vsz,vsz,N,4));
for ii=1:1:numel(center)
    tmp=ims(:,(center(ii)-vsz/2):(center(ii)+vsz/2-1),:);
    if fsign(ii)==1
        tmp=flip(tmp,2);
    end
    qds(:,:,:,ii)=tmp;
end
    