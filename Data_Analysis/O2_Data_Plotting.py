import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from pathlib import Path
import find_nearest as fn

def import_txt(name, channel, sensors):

	p = Path('/Users/Jacob/Desktop/FireStingO2_Raw/')
	file = p / name

	raw_data = np.genfromtxt(file, skip_header=14+(2*sensors-2), encoding='ISO8859', usecols = (0,1,2,3,4,5,6,7,8,9,10,11))
	data = raw_data[:,[2, 2+channel, 6+channel]]
	return data


class Experiment:

	def __init__(self, name, channel = 1, sensors = 1):
		self.name = name
		self.data = import_txt(self.name, channel, sensors)

	def plot_data(self, ax = plt, legend = True, start = 0, irradiation_start = None):

		idx = fn.find_nearest(self.data, start)[0]
		ax.plot(self.data[:,0][idx:], self.data[:,1][idx:], label = self.name, color = 'black', linewidth = 0.7)
		
		if irradiation_start is not None:
			ax.plot((irradiation_start, irradiation_start), (-0.02, 0.04), color = 'red', label = 'Irradiation Start', linewidth = 2)
		
		if legend is True:
			ax.legend()

		if ax is not plt:
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_2$ / $\mu$mol $l^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	def plot_temperature(self):

		plt.plot(self.data[:,0], self.data[:,2], label = self.name)
		plt.legend()

def main():

	exp = Experiment('200306_JS_600_4.txt', channel = 1, sensors = 1)
	#'190125_JS_526_2.txt'

	fig, ax = plt.subplots()

	exp.plot_data(ax = ax, legend = True, start = 5000., irradiation_start = 5060.)

	return fig

if __name__ == '__main__':
	main()
	plt.show()
