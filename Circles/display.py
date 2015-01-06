from pylab import *

x = linspace(-1., 1., 1001)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1, :]

sample = atleast_2d(loadtxt("sample.txt"))

figure(figsize=(8, 8))
ion()
hold(False)
  
for i in xrange(0, sample.shape[0]):
  plot([])  

  # Circle radii and positions
  width = sample[i, 5:105]
  xc = sample[i, 105:205]
  yc = sample[i, 205:305]

  # Truncate
  xc = xc[0:sample[i, 4]]
  yc = yc[0:sample[i, 4]]
  width = exp(width[0:sample[i, 4]])

  for j in xrange(0, xc.shape[0]):
    gca().add_artist(Circle((xc[j], yc[j]), width[j], alpha=0.1))
  axis([-1, 1, -1, 1])
  draw()

ioff()
show()


