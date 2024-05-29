% t: vector of frame numbers
% nfs: number of frame per section

function [cormask,fcc] = gencormask(ind,t,nfs)
fn = ceil(single(max(ind))/nfs);
v = [0:fn];
cormask = zeros(size(ind));
fcc = zeros(numel(v)-1,1);
for ii = 1:numel(v)-1
    mask = (ind >= v(ii)*nfs)&(ind < v(ii+1)*nfs);
    if sum(mask)>0
        cormask(mask) = ii;
        if ~isempty(t)
            fcc(ii) = mean(t(mask));
        else
            fcc = 0;
        end
    end
end
cormask(cormask==0) = ii;
fcc(fcc==0) = [];

end