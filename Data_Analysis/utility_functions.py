import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import io as io

def import_lamp_spectra(spectra_name):

	s = open(spectra_name).read().replace(',','.')
	data = np.genfromtxt(io.StringIO(s), skip_header = 2)

	return np.c_[data[:,0], data[:,2]], np.c_[data[:,0], data[:,3]]

def insert_image(path, x, y, zoom, ax):

	img = mpimg.imread(path)
	imagebox = OffsetImage(img, zoom = zoom)
	ab = AnnotationBbox(imagebox, (x, y), frameon = False, xycoords = ax.transAxes)
	ax.add_artist(ab)

def scientific_notation(values):

	f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
	g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.1e' % x))
	fmt = mticker.FuncFormatter(g)

	output = []

	for value in values:
		formatted = fmt(value)
		output.append(formatted)

	return output

def plot_func(data, t, plot_type, labels = ['A', 'B', 'C', 'D', 'E', 'F'], ax = plt, show_labels = False, transpose = False, markersize = 2):

	colors = ['blue', 'orange', 'green', 'red', 'black', 'grey']

	if transpose is True:
		data = data.T

	if show_labels is True:
		for counter, i in enumerate(data):
			ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = markersize, label = labels[counter])
	else:
		for counter, i in enumerate(data):
			ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = markersize)

def color_bar():

	a = np.array([[0,1]])

	plt.figure(figsize = (9, 0.5))
	img = plt.imshow(a, cmap = 'nipy_spectral')

	plt.gca().set_visible(False)
	cax = plt.axes([0.01, 0.2, 0.98, 0.6])

	cbar = plt.colorbar(orientation = 'horizontal', cax = cax)
	cbar.set_ticks([])

	plt.savefig('../Figures/Spectral_bar.pdf', transparent = True)

	plt.show()



