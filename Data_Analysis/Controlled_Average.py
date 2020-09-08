import numpy as np
import find_nearest as fn

def controlled_avg(data, values, check = True):

	raw_values = []
	idx = 0
	margin = (values[1] - values[0])/2

	if check == True:
		if values[-1] > data[:,0][-1] + margin or values[0] < data[:,0][0] - margin:
			idx_val = fn.find_nearest(values, (data[:,0][0], data[:,0][-1]))
			values = values[idx_val[0]:idx_val[1]+1]
	
	for counter, i in enumerate(values):
		raw_values.append([i,[]])
		
		while i > data[idx][0] and not abs(i - data[idx][0]) <=  margin:
			idx += 1

		tempIdx = idx

		while tempIdx < len(data) and abs(i - data[tempIdx][0]) <= margin:
			value = data[tempIdx][1]
			raw_values[counter][1].append(value)
			tempIdx += 1

	for j in raw_values:
		j[1] = sum(j[1])/len(j[1])

	avg = np.asarray(raw_values)

	return avg


