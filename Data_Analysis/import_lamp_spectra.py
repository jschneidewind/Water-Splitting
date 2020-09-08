import numpy as np
import io as io

def import_lamp_spectra(spectra_name):

	s = open(spectra_name).read().replace(',','.')
	data = np.genfromtxt(io.StringIO(s), skip_header = 2)

	return np.c_[data[:,0], data[:,2]], np.c_[data[:,0], data[:,3]]