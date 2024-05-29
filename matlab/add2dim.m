function out = add2dim(x)

sz = size(x);
if min(sz)==1
    out = zeros(1,1,numel(x));
    out(1,1,:) = x;
else
    out = zeros([1,1,sz]);
    out(1,1,:,:) = x;
end


end