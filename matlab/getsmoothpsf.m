function corrPSFhdt = getsmoothpsf(im,lambda)
b3_0t = bsarray(im,'lambda',lambda);
zhd = 1:1:b3_0t.dataSize(3);
dxxhd = 1;
[XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0t.dataSize(1),1:dxxhd:b3_0t.dataSize(2),zhd);
corrPSFhdt = interp3_0(b3_0t,XX,YY,ZZ,0);

end