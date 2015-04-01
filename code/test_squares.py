import pylab as plt
import numpy as np
import numpy.random as rng
from astropy.io import fits
from squares import *


# make real data
filename = '../data/frame-r-000094-1-0131.fits.gz'
hdulist = fits.open(filename)
realimg = hdulist[0].data
hdulist.close()
print('Done reading in %s' % filename)
realimg += 0.01* rng.normal(size=realimg.shape)


# make fake data
fakeimg = rng.normal(size=realimg.shape)
ny, nx = fakeimg.shape
ygrid, xgrid = np.mgrid[:ny, :nx]
fakeimg += 0.0 * xgrid / nx  # to add a tiny gradient as a test
startrow, startcol = 150, 350
halfrow, halfcol = np.array(fakeimg.shape)/2 + [200, 250]

tmpBBs = []
true_size = 100
fakeimg[startrow:startrow+100,startcol:startcol+true_size] += 0.2 # a patch of slightly brighter pixels.
tmpBBs.append( OddityBB([startrow,startcol], true_size) )  

true_size = 200
fakeimg[halfrow:halfrow+true_size,startcol:startcol+true_size] *= 1.2 # a patch of slightly MORE variable pixels.
tmpBBs.append( OddityBB([halfrow,startcol], true_size) )  

true_size = 200
fakeimg[halfrow:halfrow+true_size,halfcol:halfcol+true_size] *= 0.7 # a patch of slightly less variable pixels.
tmpBBs.append( OddityBB([halfrow,halfcol], true_size) )  

true_size = 250
fakeimg[startrow:startrow+true_size,halfcol:halfcol+true_size] -= 0.2 # a patch of slightly darker pixels.
tmpBBs.append( OddityBB([startrow,halfcol], true_size) )  

show_data_and_model(fakeimg, 'truth', tmpBBs)
print('Done inventing fake data')

img = fakeimg # we decide which data we're working with.
show_data_and_model(img, 'raw')

# ------------------------------------------------------------------

M = 32 # number of histogram bins


cbb = intensity_to_cumulated_block(img, M)

# initialising MCMC with some BBs
halfrow, halfcol = np.array(img.shape)/2
BBs = []
initsize = min(halfrow,halfcol) # deliberately so the BBs are huge but each enclose one true source.
BBs.append( OddityBB([0,0], initsize) )  
BBs.append( OddityBB([0,initsize], initsize) )  
BBs.append( OddityBB([initsize,0], initsize) )  
BBs.append( OddityBB([initsize,initsize], initsize) )  

#show_data_and_model(intensity_to_boolean_block(img, M)[-1], 'raw_top_bin')
show_data_and_model(img, 'start', BBs)

alphas_S = np.ones(M, dtype=int) 
alphas_B = 1 * alphas_S


indices = range(len(BBs))
logP = MH_step(cbb, 0, BBs, alphas_S, alphas_B)
T = 0
for t in range(7):
    T_inner = 2 ** (2*t)
    for tt in range(T_inner):
        for i in indices:
            logP = MH_step(cbb, i, BBs, alphas_S, alphas_B, temperature=1., logP=logP)
    T += T_inner
    print '%5d' % (T), 
    for i in indices:
        print BBs[i],
    print '%.1f' % (logP)

show_data_and_model(img, 'end', BBs)

top_bin_img = intensity_to_boolean_block(img, M)[-1]
show_data_and_model(top_bin_img, 'top_bin', BBs)

# show the histograms of each BB, and the one for the background



plt.clf()
sources_counts = 0
max_c = 0
for i, BB in enumerate(BBs):
    this_counts = BB.get_counts(cbb)
    this_logL_contribution = get_logL_patch(this_counts, alphas_S)
    sources_counts += this_counts
    max_c = max(max_c,  np.max(this_counts))
    plt.subplot(1, len(BBs)+1, i+1)
    plt.plot(this_counts, 'ok', this_counts, '-k')
    plt.gca().set_ylim([0,50000])
    if i>0: plt.gca().set_yticks([])

for i, BB in enumerate(BBs):
    plt.subplot(1, len(BBs)+1, i+1)
    plt.gca().set_ylim([0,1.1*max_c])
    plt.title('source %i' %(i) )

background_counts = cbb[:, -1, -1].ravel() - sources_counts
BG_logL_contribution = get_logL_patch(background_counts, alphas_B)
plt.subplot(1, len(BBs)+1, len(BBs)+1)
plt.plot(background_counts, 'or', background_counts, '-r')
plt.gca().set_ylim([0,np.max(background_counts)])
plt.title('background')
plt.savefig('histograms.png',dpi=150)

plt.clf()
plt.hist(np.ravel(img))
plt.show()
