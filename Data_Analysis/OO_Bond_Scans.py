import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

class OO_Bond_Scan:

	def __init__(self, name, de_scan = False):
		self.name = name

		if de_scan is False:
			self.import_data()
		else:
			self.import_de_data()

	def import_de_data(self):
		self.df = pd.read_excel('../Computational_Data/DE_Scan/%s' % self.name, dtype = object)
		self.dict = self.df.to_dict()

	def import_data(self):
		self.df = pd.read_excel('../Computational_Data/Misc/OO_Bond_Scans/OO_Bond_Scans.xlsx', dtype = object, sheet_name = self.name)
		self.dict = self.df.to_dict()

	def plot_image(self, x_lim, y_lim, ax, path = '../Computational_Data/Images/', i = 0, xycoords = 'data'):

		image = self.dict['Image'][i]

		if image != 'No':

			vert_off = self.dict['Vertical_Offset'][i]
			horiz_off = self.dict['Horizontal_Offset'][i]

			img = mpimg.imread(path + image)
			imagebox = OffsetImage(img, zoom = self.dict['Zoom'][i])

			ab = AnnotationBbox(imagebox, (x_lim - horiz_off, y_lim + vert_off), frameon = False, xycoords = xycoords)
			ax.add_artist(ab)

	def plot_data(self, ax = plt, label = None, i = 0, bond_length = 'Bond Length', energy = 'Energy (kcal/mol)', xlabel = r'O-O Distance / $\AA$', x_lim_upper = 2.95, x_lim_lower = 1.3, image_path = '../Computational_Data/Images/', plot_image = True, markersize = 10., plot_grid = True):

		if label is None:
			label = self.dict['Name'][i]

		ax.plot(self.df[bond_length], self.df[energy], 'o-', linewidth = 2., markersize = markersize, label = label, zorder = 10)

		if ax != plt:

			ax.set_xlabel(xlabel)
			ax.set_ylabel(r'Relative Energy / kcal mol$^{-1}$')
			ax.legend()

			ax.set_xlim(x_lim_upper, x_lim_lower)
			y_lim = ax.get_ylim()[1]

			if plot_image is True:
				self.plot_image(x_lim_upper, y_lim, ax, path = image_path, i = i)

			if plot_grid is True:
				ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

def main():

	gme_cis = OO_Bond_Scan('GMe-Cis-Mono-D-0')
	gme = OO_Bond_Scan('GMe-Mono-D-0')
	three_S = OO_Bond_Scan('3-Full-Disp')
	three_T = OO_Bond_Scan('3-Full-T')
	GAlt_S = OO_Bond_Scan('GAlt-Mono-S-0')
	GAlt_T = OO_Bond_Scan('GAlt-Mono-T-0')

	fig = plt.figure(figsize = (11, 7))

	gs = gridspec.GridSpec(2, 4)
	gs.update(wspace=0.35)
	ax0 = plt.subplot(gs[0, :2], )
	ax1 = plt.subplot(gs[0, 2:])
	ax2 = plt.subplot(gs[1, 1:3])

	fig.subplots_adjust(hspace = 0.35, top = 0.92, bottom = 0.07, left = 0.1, right = 0.9)

	gme_cis.plot_data(ax = ax1)
	gme.plot_data(ax = ax1)
	three_S.plot_data(ax = ax0)
	three_T.plot_data(ax = ax0)
	GAlt_S.plot_data(ax = ax2)
	GAlt_T.plot_data(ax = ax2)

	ax0.text(-0.07, 1.05, 'A', transform=ax0.transAxes, size = 20, weight='bold')
	ax1.text(-0.07, 1.05, 'B', transform=ax1.transAxes, size = 20, weight='bold')
	ax2.text(-0.07, 1.05, 'C', transform=ax2.transAxes, size = 20, weight='bold')

	return fig

def secondary():

	de_scan = OO_Bond_Scan('DE_Scan_Energies.xlsx', de_scan = True)
	
	fig, ax = plt.subplots(figsize = (5, 7))
	fig.subplots_adjust(bottom = 0.4, top = 0.6)

	points = 7.
	x_lim_upper = 2.65
	x_lim_lower = 6.
	image_pos = 0.5
	arrow_cut_off = 0.5

	de_scan.plot_data(ax = ax, bond_length = 'Bond Length (A)', energy = 'Gas Phase Energy (kcal/mol)', x_lim_upper = x_lim_upper, x_lim_lower = x_lim_lower, 
					  xlabel = r'Ru-O Distance / $\AA$', plot_image = False, label = 'Gas Phase Energy', markersize = 7, plot_grid = True)

	de_scan.plot_data(ax = ax, bond_length = 'Bond Length (A)', energy = 'SMD Energy (kcal/mol)', x_lim_upper = x_lim_upper, x_lim_lower = x_lim_lower, 
					  xlabel = r'Ru-O Distance / $\AA$', plot_image = False, label = 'SMD Energy', markersize = 7, plot_grid = True)

	for i in de_scan.dict['Image']:
		de_scan.plot_image(i/points, image_pos, ax, path = '../Computational_Data/DE_Scan/J-Trans-Mono-H2O-T0-RuO-Scan-GPO/', i = i, xycoords = ax.transAxes)

		bond_length = de_scan.dict['Bond Length (A)'][i]
		line_position = (bond_length - x_lim_upper)/(x_lim_lower - x_lim_upper)
		vert_off = de_scan.dict['Vertical_Offset'][i]

		dotted_line = ax.annotate('', xy = (line_position, 0), xytext = (line_position, 1), arrowprops={'arrowstyle': '-', 'ls': 'dashed', 'color': 'grey'}, xycoords = ax.transAxes, textcoords = ax.transAxes, zorder = 1)

		if i == 3:
			arrow = ax.annotate('', xy = (line_position, image_pos + vert_off + arrow_cut_off + 0.2), xytext = (line_position, -0.05), arrowprops={'arrowstyle': '->', 'ls': 'dashed', 'color': 'grey'}, xycoords = ax.transAxes, textcoords = ax.transAxes, zorder = 0)

		else:

			if vert_off > 0:
				arrow = ax.annotate('', xy = (line_position, image_pos + vert_off - arrow_cut_off), xytext = (line_position, 1.05), arrowprops={'arrowstyle': '->', 'ls': 'dashed', 'color': 'grey'}, xycoords = ax.transAxes, textcoords = ax.transAxes, zorder = 0)

			else:
				arrow = ax.annotate('', xy = (line_position, image_pos + vert_off + arrow_cut_off), xytext = (line_position, -0.05), arrowprops={'arrowstyle': '->', 'ls': 'dashed', 'color': 'grey'}, xycoords = ax.transAxes, textcoords = ax.transAxes, zorder = 0)

	return fig

if __name__ == '__main__':
	#main()
	secondary()
	plt.show()
