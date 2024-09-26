#%% -*- coding: utf-8 -*-
"""
3D drift estimation example
Units for X,Y,Z are pixels, pixels, and microns resp.
"""
import numpy as np
import matplotlib.pyplot as plt

from dme.dme import dme_estimate
from dme.rcc import rcc3D
from dme.native_api import NativeAPI
import glob
import os
import scipy.io as sio
import h5py as h5
# Need to have CUDA installed
use_cuda=True
#%%
folder = r'E:\EMBL files\data for PSF learning\02-18-2022 nup96 mMaple\psf insitu fit/'

filename = 'cell3_nup96_mMaple__560_00001_00000_loc_z_smap.mat'

resfile = folder+filename
F = h5.File(resfile,'r')

res = F['res']
g = np.array([129,129,1000])
localizations = np.swapaxes(np.vstack([res['x'],res['y'],res['z']]),0,1)/g
crlb = np.swapaxes(np.vstack([res['stdx'],res['stdy'],res['stdx']]),0,1)/g # standard deviation
framenum = np.squeeze(np.int64(res['frames']))
fov_width = 200


#%% Simulate an SMLM dataset in 3D with blinking molecules

print(f"Total localizations: {len(localizations)}")
    
estimated_drift = dme_estimate(localizations, framenum, 
             crlb, 
             framesperbin = 100,  # note that small frames per bin use many more iterations
             imgshape=[fov_width, fov_width], 
             coarseFramesPerBin=1000,
             coarseSigma=[0.2,0.2,0.2],  # run a coarse drift correction with large Z sigma
             rccZoom=4,
             estimatePrecision=False,
             useCuda=use_cuda,
             useDebugLibrary=False)
#%%
estimated_drift_rcc = rcc3D(localizations, framenum, timebins=10, zoom=4)

#rmsd = np.sqrt(np.mean((estimated_drift-drift_trace)**2, 0))
#print(f"RMSD of drift estimate compared to true drift: {rmsd}")

fig,ax=plt.subplots(3, figsize=(7,6))
for i in range(3):
    #ax[i].plot(drift_trace[:,i],label='True drift')
    ax[i].plot(estimated_drift[:,i]+0.2,label='Estimated drift (DME)')
    ax[i].plot(estimated_drift_rcc[:,i]-0.2,label='Estimated drift (RCC)')
    ax[i].set_title(['x', 'y', 'z'][i])

    unit = ['px', 'px', 'um'][i]
    ax[i].set_ylabel(f'Drift [{unit}]')
ax[0].legend()
plt.tight_layout()
plt.show()

# %%
res1 = {'shift_dme':estimated_drift,'shift_rcc':estimated_drift_rcc}

matfile = os.path.splitext(resfile)[0]+'_dmeshift.mat'
sio.savemat(matfile,res1)
# %%
