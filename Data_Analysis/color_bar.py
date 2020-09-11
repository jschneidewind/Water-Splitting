import matplotlib.pyplot as plt
import numpy as np

a = np.array([[0,1]])

plt.figure(figsize = (9, 0.5))
img = plt.imshow(a, cmap = 'nipy_spectral')

plt.gca().set_visible(False)
cax = plt.axes([0.01, 0.2, 0.98, 0.6])

cbar = plt.colorbar(orientation = 'horizontal', cax = cax)
cbar.set_ticks([])

plt.savefig('../Figures/Spectral_bar.pdf', transparent = True)

plt.show()