import numpy as np
import numpy.random as rng

class Squares:
  # Image dimensions
  x_min = -1.
  x_max =  1.
  y_min = -1.
  y_max =  1.
  x_range = x_max - x_min
  y_range = y_max - y_min
  scale = np.sqrt(x_range*y_range)

  # Maximum number of squares
  N_max = 1000

  def __init__(self):
    """
    Create model with zero squares
    """
    self.xc = np.empty(Squares.N_max)
    self.yc = np.empty(Squares.N_max)
    self.widths = np.empty(Squares.N_max)
    self.N = 0

  def add(self):
    """
    Add a square to the model
    Return number of collisions created
    """
    if self.N == Squares.N_max:
      return 0

    self.xc[self.N] = Squares.x_min + Squares.x_range*rng.rand()
    self.yc[self.N] = Squares.y_min + Squares.y_range*rng.rand()
    self.widths[self.N] = 1E-2*Squares.scale

    diffx = self.xc[self.N] - self.xc[0:self.N]
    diffy = self.yc[self.N] - self.yc[0:self.N]
    collision = np.logical_and(\
      diffx < self.widths[self.N] + self.widths[0:self.N], \
      diffy < self.widths[self.N] + self.widths[0:self.N])

    self.N += 1
    return np.sum(collision)

  def remove(self):
    """
    Remove a square from the model
    """
    if self.N == 0:
      return

    which = rng.randint(self.N)
    self.xc[which:(self.N-1)] = self.xc[(which+1):self.N]
    self.yc[which:(self.N-1)] = self.yc[(which+1):self.N]
    self.widths[which:(self.N-1)] = self.widths[(which+1):self.N]
    self.N -= 1

  def proposal(self):
    if rng.rand() < 0.5:
      self.add()
    else:
      self.remove()

if __name__ == '__main__':
  from pylab import *
  from copy import deepcopy

  x = linspace(-1., 1., 1001)
  y = x.copy()
  [x, y] = meshgrid(x, y)
  y = y[::-1, :]

  figure(figsize=(8, 8))
  ion()
  hold(False)

  s = Squares()

  for i in xrange(0, 1000):
    plot([])
    s.proposal()

    # Circle radii and positions
    width = s.widths[0:s.N]
    xc = s.xc[0:s.N]
    yc = s.yc[0:s.N]

    for j in xrange(0, xc.shape[0]):
      gca().add_artist(Rectangle((xc[j] - width[j], yc[j] - width[j]), 2*width[j], 2*width[j], alpha=0.1))
    axis([-1, 1, -1, 1])
    title(i+1)
    draw()

  ioff()
  show()


