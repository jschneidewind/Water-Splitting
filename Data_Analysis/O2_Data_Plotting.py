import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from pathlib import Path

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

	def plot_data(self):

		plt.plot(self.data[:,0], self.data[:,1], label = self.name)
		plt.legend()

	def plot_temperature(self):

		plt.plot(self.data[:,0], self.data[:,2], label = self.name)
		plt.legend()



exp = Experiment('200306_JS_600_4.txt', channel = 1, sensors = 1)
#'190125_JS_526_2.txt'



exp.plot_data()
exp.plot_temperature()

plt.show()
