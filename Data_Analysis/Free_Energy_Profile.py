import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def get_name(name_label):

	name_split = str.split(name_label, ';')
	name_string = ''

	for counter_name, j in enumerate(name_split):
		name_string += j
		if len(name_split) > 1 and counter_name < len(name_split) - 1:
			name_string += '\n'

	return name_string

class Free_Energy_Profile:

	def __init__(self, name):
		self.name = name
		self.import_data()
		self.final_figure = self.plot_profile()

	def import_data(self):

		df = pd.read_excel(self.name, index_col = 0, dtype = object)

		self.energies = df['Relative Energy (kcal/mol)'].to_numpy()
		self.data = df.to_dict('index')		

	def plot_plus_sign(self, int_type, target, pos, energy, plus_sign_offset, ax):

		if int_type == target:
			plus = ax.annotate('+', xy = (pos, energy + plus_sign_offset), fontsize = 15, ha = 'center')

	def plot_seperator(self, seperator, pos, y_lim_low, seperator_bottom, energy, seperator_energy, seperator_top_middle, y_lim_high, seperator_top, ax):

		if seperator != 'No':
			if seperator == 'bottom':
				seperation = ax.annotate('', xy = (pos, y_lim_low + seperator_bottom), xytext = (pos, energy - seperator_energy), arrowprops={'arrowstyle': '-', 'ls': 'dashed', 'color': 'grey'})
			else:
				seperation = ax.annotate('', xy = (pos, energy + seperator_energy + seperator_top_middle), xytext = (pos, y_lim_high + seperator_top), arrowprops={'arrowstyle': '-', 'ls': 'dashed', 'color': 'grey'})

	def plot_connection_excitation(self, counter, int_type, prev_int_type, prev_pos, width, prev_energy, pos, energy, arrow_space, ax, color_used):

		if counter > 0 and int_type != prev_int_type:

			if int_type == 'B' or int_type == 'F':
				connection = ax.annotate('', xy =(prev_pos + width, prev_energy), xytext=(pos - width, energy), arrowprops={'arrowstyle': '-', 'ls': 'dotted', 'color': 'grey'})

			else:
				connection = ax.annotate('', xy =(prev_pos + width, prev_energy), xytext=(pos - width, energy), arrowprops={'arrowstyle': '-', 'ls': 'dashed', 'color': 'black'})

		elif int_type == prev_int_type:
			if int_type == 'A':
				wavelength = ax.annotate('h$\\nu$\n(320-\n400 nm)', xy = (pos + self.wavelength_offset, (energy + prev_energy)/2.), fontsize = 8, ha = 'center', va = 'center')
				excitation = ax.annotate('', xy =(pos, energy - arrow_space), xytext=(pos, prev_energy + arrow_space), arrowprops={'arrowstyle': '->', 'color': 'darkblue'})
			if int_type == 'B':
				wavelength = ax.annotate('h$\\nu$\n(450-\n600 nm)', xy = (pos + self.wavelength_offset, (energy + prev_energy)/2.), fontsize = 8, ha = 'center', va = 'center')
				excitation = ax.annotate('', xy =(pos, energy - arrow_space), xytext=(pos, prev_energy + arrow_space), arrowprops={'arrowstyle': '->', 'color': 'darkred'})

	def plot_line_and_label(self, pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = 0, energy_label = True, intermediate_label = True, approx = False, color_text = 'black'):

		line = ax.annotate('', xy =(pos - width, energy), xytext=(pos + width, energy), arrowprops={'width': 2.5, 'color': color_used, 'headwidth': 1, 'headlength': 0.005})
		
		if intermediate_label == True:
			label = ax.annotate(r'%s' %name_label, xy =(pos, energy + text_offset), ha = 'center')

		if energy_label == True:

			if approx == True:
				energy_label = ax.annotate('~ {:.0f}'.format(energy), xy = (pos, energy - energy_offset), ha = 'center')

			else:
				energy_label = ax.annotate('{:.1f}'.format(energy), xy = (pos, energy - energy_offset), ha = 'center', color = color_text)

	def plot_intermediate(self, i, prev_energy, int_type, prev_int_type, pos, width, energy, color_used, name_label, text_offset, energy_offset, offset_multi, ax):

		if int_type == 'A' or int_type == 'B':

			if int_type == prev_int_type:
				self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = energy_offset, approx = True)
			else:	
				self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = energy_offset)

		elif int_type != prev_int_type and int_type != 'BC' and int_type != 'D':

			if i == 'DE_Scan_5':
				self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = energy_offset, approx = True)

			else:
				self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = energy_offset)

		elif int_type == 'BC':
			ax.plot(pos, energy, '.', markersize = 35, color = color_used)
			self.plot_line_and_label(pos, width, energy, color_used, name_label, offset_multi * text_offset, ax, energy_label = False)

		elif i == 'D':
			self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_label = False)

		elif i == 'D_Mono':
			self.plot_line_and_label(pos, width, energy, color_used, name_label, text_offset, ax, energy_offset = energy_offset, intermediate_label = False, energy_label = False)
			energy_label = ax.annotate('{:.1f}'.format(prev_energy), xy = (pos, energy - energy_offset), ha = 'center', color = 'darkblue')
			energy_label = ax.annotate('{:.1f}'.format(energy), xy = (pos, energy - 2 * energy_offset), ha = 'center', color = 'darkgreen')

	def plot_image(self, i, pos, energy, image_offset, ax):

		image = self.data[i]['Image']

		if image != 'No':

			vert_off = self.data[i]['Vertical_Offset']
			horiz_off = self.data[i]['Horizontal_Offset']

			img = mpimg.imread('../Computational_Data/Images/%s' % image)
			imagebox = OffsetImage(img, zoom = self.data[i]['Zoom'])

			if self.data[i]['Top_Bottom'] == 'top':
				ab = AnnotationBbox(imagebox, (pos + horiz_off, energy + image_offset + vert_off), frameon = False)
				ax.add_artist(ab)

			else:
				ab = AnnotationBbox(imagebox, (pos + horiz_off, energy - image_offset + vert_off), frameon = False)
				ax.add_artist(ab)

	def plot_oxidation_state(self, i, ox_offset, ox_sphere_zoom, pos, ax):

		ox_state = self.data[i]['Oxidation_State']
		ox_offsets = {0: ox_offset, 1: -ox_offset}

		if ox_state != 'No':

			ox_split = str.split(ox_state, ' ')

			for counter, i in enumerate(ox_split):

				sphere = mpimg.imread('../Computational_Data/Images/%s' % i)
				spherebox = OffsetImage(sphere, zoom = ox_sphere_zoom)

				if len(ox_split) == 1:
					sphere_ab = AnnotationBbox(spherebox, (pos, 0), frameon = False)
					ax.add_artist(sphere_ab)

				else:
					sphere_ab = AnnotationBbox(spherebox, (pos + ox_offsets[counter], 0), frameon = False)
					ax.add_artist(sphere_ab)

	def plot_profile(self):

		fig, ax = plt.subplots(figsize = (11, 6.5))
		fig.subplots_adjust(wspace = 0, hspace = 0 )
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.set_xticks([]) 
		ax.set_xlabel('Reaction Coordinate')
		ax.set_ylabel(r'$\Delta$G / kcal mol$^{-1}$')
		color_model_chemistry = {'Dimer': 'darkblue', 'Monomer (CASSCF)': 'darkred', 'Monomer': 'darkgreen'}
 
		x_lim_low = -1.
		x_lim_high = 9
		y_lim_low = -10.
		y_lim_high = 139.
		position = {'A': -0.5, 'B': 1, 'BC': 2, 'C': 3, 'CD': 4, 'D': 5, 'DE': 6, 'E': 7, 'F': 8.5}
		
		width = 0.3
		arrow_space = 7.
		text_offset = 3.
		energy_offset = 5.
		image_offset = 20.
		offset_multi = 2.5
		seperator_top = 10.
		seperator_top_middle = 10.
		seperator_energy = 10.
		seperator_bottom = 40.
		ox_offset = 0.35
		ox_sphere_zoom = 0.023
		plus_sign_offset = -25.
		self.wavelength_offset = 0.3

		ax.set_xlim(x_lim_low, x_lim_high)
		ax.set_ylim(y_lim_low, y_lim_high)
		ax.plot((x_lim_high), (y_lim_low), ls="", marker=">", ms=5, color="k", clip_on=False)   
		ax.plot((x_lim_low), (y_lim_high), ls="", marker="^", ms=5, color="k", clip_on=False)
           
		index_list = list(self.data)

		for counter, i in enumerate(self.data):

			energy = self.data[i]['Relative Energy (kcal/mol)']
			model = self.data[i]['Model_Chemistry']
			name_label = self.data[i]['Name']
			seperator = self.data[i]['Line']
			color_used = color_model_chemistry[model]

			name_label = get_name(name_label)
			int_type = str.split(i, '_')[0]
			pos = position[int_type]

			ax.plot(pos, energy, '-', color = color_used, label = model, markersize = 10, mew = 20, linewidth = 4)

			prev_int = index_list[counter-1]
			prev_int_type = str.split(prev_int, '_')[0]
			prev_pos = position[prev_int_type]
			prev_energy = self.data[prev_int]['Relative Energy (kcal/mol)']

			self.plot_plus_sign(int_type, 'E', pos, energy, plus_sign_offset, ax)

			self.plot_seperator(seperator, pos, y_lim_low, seperator_bottom, energy, seperator_energy, seperator_top_middle, y_lim_high, seperator_top, ax)

			self.plot_connection_excitation(counter, int_type, prev_int_type, prev_pos, width, prev_energy, pos, energy, arrow_space, ax, color_used)

			self.plot_intermediate(i, prev_energy, int_type, prev_int_type, pos, width, energy, color_used, name_label, text_offset, energy_offset, offset_multi, ax)

			self.plot_image(i, pos, energy, image_offset, ax)

			self.plot_oxidation_state(i, ox_offset, ox_sphere_zoom, pos, ax)

		handles, labels = ax.get_legend_handles_labels()
		by_label = dict(zip(labels, handles))
		ax.legend(by_label.values(), by_label.keys(), loc = (0.78, 0.17), title = r'$\bf{Model}$ $\bf{Chemistry}$')

		#fig.savefig('Test_save_new.pdf', dpi = 500, transparent = True)

		return fig

class Conical_Intersections:

	def __init__(self, name, pcm = False):
		self.pcm = pcm
		self.name = name
		self.import_data(pcm = pcm)
		self.final_figure = self.plot_CIs()

	def import_data(self, pcm):

		if pcm is False: 
			df = pd.read_excel(self.name, sheet_name = 0, index_col = 0, dtype = object)
		else:
			df = pd.read_excel(self.name, sheet_name = 1, index_col = 0, dtype = object)

		ints = pd.read_excel(self.name, sheet_name = 2, index_col = 0, dtype = object)

		self.data = df.to_dict()
		self.data_numpy = df.to_numpy()
		self.intermediates = ints.to_dict('index')

	def insert_image(self, roots, image, zoom, counter, image_offset, image_offset_y, horiz_off, vert_off, ax):

		img = mpimg.imread('../Computational_Data/Images/%s' % image)
		imagebox = OffsetImage(img, zoom = zoom)		
		fraction_ax = float(counter)/float(roots)
		ab = AnnotationBbox(imagebox, (fraction_ax - image_offset + horiz_off, image_offset_y + vert_off), frameon = False, xycoords = ax.transAxes)
		ax.add_artist(ab)

	def plot_CIs(self):

		fig, ax = plt.subplots(figsize = (5, 6.5))

		x_labels = []
		fig.subplots_adjust(bottom=0.5)
		fig.subplots_adjust(top=0.95)

		positions = np.array([0, 1, 2, 3])

		x_lim_low = -0.1
		x_lim_high = 3.1
		y_lim_low = -5.

		if self.pcm is False:
			y_lim_high = 105.
		else:
			y_lim_high = 110.

		roots = 3
		idx = 5 + roots
		path_offset = 2.9
		image_offset = 0.1
		image_offset_y = -0.4
		sep_offset = -0.13
		sep_end = -0.68
		vert_sep_pos = -0.63

		ax.set_xlim(x_lim_low, x_lim_high)
		ax.set_ylim(y_lim_low, y_lim_high)
		ax.set_ylabel(r'Relative Energy / kcal mol$^{-1}$')

		roots_names = ('VDZ kcal Root 1', 'VDZ kcal Root 2', 'VDZ kcal Root 3')
		data = self.data_numpy[5:idx]

		path = np.array([data[2][0], data[2][1], data[1][2], data[0][3]])

		for counter, i in enumerate(self.intermediates):

			name_label = self.intermediates[i]['Full_Name']
			name_label = get_name(name_label)
			x_labels.append(name_label)

			dotted_line = ax.annotate('', xy = (counter, y_lim_low), xytext = (counter, y_lim_high), arrowprops={'arrowstyle': '-', 'ls': 'dashed', 'color': 'grey'})

			image = self.intermediates[i]['Image']
			vert_off = self.intermediates[i]['Vertical_Offset']
			horiz_off = self.intermediates[i]['Horizontal_Offset']

			self.insert_image(roots, image, self.intermediates[i]['Zoom'], counter, image_offset, image_offset_y, horiz_off, vert_off, ax)			

			sep_position = (counter + abs(x_lim_low))/(abs(x_lim_low) + x_lim_high)

			if counter == 1 or counter == 3:
				vert_seperator = ax.annotate('', xy = (sep_position, sep_end), xytext = (sep_position, sep_offset), xycoords = ax.transAxes, textcoords = ax.transAxes, arrowprops={'arrowstyle': '->', 'ls': 'dashed', 'color': 'grey'}, zorder = 1)

		for counter, i in enumerate(data[::-1]):
			ax.plot(positions, i, 'o-', label = r'$\bf{D_{%s}}$' % (roots - counter - 1), linewidth = 2.5, markersize = 10, zorder = 5)
		
		ax.plot(positions, path + path_offset, '--', color = 'black', label = 'Reaction Path', zorder = 10)

		plt.xticks(positions, x_labels)
		ax.legend()

		return fig

		#fig.savefig('Test_CI.pdf', dpi = 500, transparent = True)

if __name__ == '__main__':

	profile = Free_Energy_Profile('../Computational_Data/Overview_Energies/DFT_Energies.xlsx')
	#conical = Conical_Intersections('../Computational_Data/Overview_Energies/CASSCF_Energies.xlsx', pcm = False)

	plt.show()