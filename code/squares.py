import pylab as plt
import numpy as np
import numpy.random as rng
from astropy.io import fits


def show_data_and_model(BG_image, title='dunno', BoundingBoxes = []):
    plt.clf()
    plt.imshow(BG_image, interpolation='nearest', origin='lower')
    for BB in BoundingBoxes:
        plt.gca().add_artist(plt.Rectangle((BB.get_pos()[1], BB.get_pos()[0]), BB.get_width(), BB.get_width(), alpha=0.9, facecolor='None', edgecolor='blue'))
    plt.title(title)
    plt.savefig(title+'.png',dpi=150)
    print('Done writing %s' % (title+'.png'))

class OddityBB:
    def __init__(self, pos, width):
        """
        Warning: this sets pos and width to be ints; truncation.
        """
        self.set_pos(pos)
        self.set_width(width)
    
    def set_pos(self, pos):
        assert len(pos) == 2
        self.pos = np.array(pos, dtype=int)
        
    def set_width(self, width):
        self.width = int(width)
        assert self.width > 0

    def __str__(self):
        return "Oddity BB, pos %4d %4d, width %4d" % (self.get_pos()[0], self.get_pos()[1], self.get_width())
    
    def __eq__(self, other):
        if np.all(self.get_pos() == other.get_pos()) and self.get_width() == other.get_width():
            return True
        return False

    def is_overlapping(self, other):
        dy, dx = self.pos - other.pos
        if dy > 0.:
            ywidth = other.width
        else:
            ywidth = self.width
        if dx > 0.:
            xwidth = other.width
        else:
            xwidth = self.width
        return np.abs(dy) < ywidth and np.abs(dx) < xwidth
    
    def _get_corners(self):
        """
        We're defining what we mean by corners here! 
        These are the 'corners' relevant to the output of intensity_to_cumulated_block.
        """
        return (self.pos,
                self.pos + [self.width, 0],
                self.pos + [self.width, self.width],
                self.pos + [0, self.width])
    
    def get_pos(self):
        return self.pos  
    
    def get_width(self):
        return self.width
    
    def _get_cumulated_block_elt(self, pos, cumulated_block):
        M, nyplus1, nxplus1 = cumulated_block.shape
        cpos = np.clip(pos, [0, 0], [nyplus1-1, nxplus1-1])
        return cumulated_block[:, cpos[0], cpos[1]]
    
    def get_counts(self, cumulated_block):
        corners = self._get_corners()
        return (self._get_cumulated_block_elt(corners[0], cumulated_block)
                - self._get_cumulated_block_elt(corners[1], cumulated_block)
                + self._get_cumulated_block_elt(corners[2], cumulated_block)
                - self._get_cumulated_block_elt(corners[3], cumulated_block))



def intensity_to_m(img, M):
    ny, nx = img.shape
    N = ny*nx
    mimgvec = np.zeros(N,dtype=int)
    mimgvec[np.argsort(img.reshape(N))] = np.floor((float(M) / float(N)) * np.arange(N)).astype(int)
    return mimgvec.reshape((ny, nx))


def intensity_to_boolean_block(img, M):
    mimg = intensity_to_m(img, M)
    ny, nx = mimg.shape
    N = ny * nx
    bb = np.zeros((M, ny+1, nx+1),dtype=bool)
    for m in range(M):
        bb[m, 1:, 1:] = (mimg == m)
    return bb


def get_logL_patch(counts, alphas):
    """
    Takes a vector of counts (across bins), and vector of alpha hyperparameters (ditto).
    Returns the log likelihood of those counts under the Dirichlet-multinomial distribution with
    those hyperparameters.
    """
    # initialize internal (static) variables
    try:
        get_logL_patch.Cmax
    except:
        get_logL_patch.Cmax = 1
        get_logL_patch.gammaln = np.array([np.Inf])
    N, A = np.sum(counts), np.sum(alphas)
    # calculate more lookup table if necessary
    if N + A > get_logL_patch.Cmax:
        print "get_logL_patch(): Calculating expensive lookup shit"
        get_logL_patch.Cmax = 2 * (N + A)
        get_logL_patch.gammaln = np.append(np.array([np.Inf, 0.]), np.cumsum(np.log(np.arange(get_logL_patch.Cmax - 2) + 1)))
    # now the actual LF
    return (  get_logL_patch.gammaln[A] 
            - get_logL_patch.gammaln[N + A] 
            + np.sum(get_logL_patch.gammaln[counts + alphas]) 
            - np.sum(get_logL_patch.gammaln[alphas]) 
            )


def intensity_to_cumulated_block(img, M):
    return np.cumsum(np.cumsum(intensity_to_boolean_block(img, M), axis=2), axis=1)


def get_logL(cumulated_block, BBs, alphas_S, alphas_B):
    """
    ## inputs
    - `BBs` is a list of `OddityBB` objects.
    """
    logL = 0.
    source_counts = 0
    for BB in BBs:
        this_counts = BB.get_counts(cumulated_block)
        logL += get_logL_patch(this_counts, alphas_S)
        source_counts += this_counts
    background_counts = cumulated_block[:, -1, -1].ravel() - source_counts
    logL += get_logL_patch(background_counts, alphas_B)
    return logL


def get_logPrior(BBs):
    logP = 0.
    if len(BBs) == 0:
        return logP
    for i in range(len(BBs)):
        logP += -2. # MAGIC NUMBER, the cost of adding a new source.
        if BBs[i].get_width() > 500: # MAGIC NUMBER, the biggest "source" we believe in.
            logP += -1. * (BBs[i].width - 500)
        if BBs[i].get_width() < 9: # MAGIC NUMBER, the minimum source size we find interesting
            #print 'get_logPrior() : bounding box %s is too narrow' %(BBs[i])
            return -np.Inf
    for i in range(len(BBs) - 1):
        for j in range(i + 1, len(BBs)):
            if BBs[i].is_overlapping(BBs[j]):
                #print 'get_logPrior() : bounding boxes %s and %s overlap' %(BBs[i], BBs[j])
                return -np.Inf
    return logP



def get_logPosterior(cbb, BBs, alphas_S, alphas_B):
    logP = get_logPrior(BBs)
    if np.isinf(logP):
        #print 'get_logPosterior() : infinite logP' 
        return -np.Inf
    return logP + get_logL(cbb, BBs, alphas_S, alphas_B)


def MH_step(cbb, i, BBs, alphas_S, alphas_B, temperature = 1., logP = None):
    """
    Metropolis Hastings move on BBs[i]; return the next
    Also change BBs in place, so FUCKING BE CAREFUL.
    """
    if logP == None:
        logP = get_logPosterior(cbb, BBs, alphas_S, alphas_B)
    pos = BBs[i].get_pos()
    width = BBs[i].get_width()
    jump_size = 2 ** rng.randint(np.log2(np.min(cbb[0].shape)/8))
    # 1, 2, 4, and so on... 
    new_pos = pos + jump_size * rng.randint(-1,2,size=2) 
    # but new_pos should stay in bounds 
    new_pos = np.maximum(np.minimum(new_pos, cbb[0].shape), [0,0]) 
    new_width =  min(max(width + jump_size * rng.randint(-1,2),1),min(cbb[0].shape)) # keep positive, girls
    BBs[i].set_pos(new_pos)
    BBs[i].set_width(new_width)
    new_logP = get_logPosterior(cbb, BBs, alphas_S, alphas_B)
    # M-H step: 
    # We should accept proposal with probability P_accept = min(1, newP / P).
    # eg. make a random uniform in [0,1], and accept if it is < exp(new_logP - logP)).
    # Or take log of the random and accept if that is < (new_logP - logP).
    # But we've ALREADY made the change, so we have to UNDO this if we reject.
    if not(np.log(rng.uniform()) < (new_logP - logP)/temperature): 
        # Reject the proposal, so go back to the old values.
        BBs[i].set_pos(pos)
        BBs[i].set_width(width)
        return logP
    return new_logP

