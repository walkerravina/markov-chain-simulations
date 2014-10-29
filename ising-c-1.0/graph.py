import matplotlib.pyplot as plt
import sys

def main():
	make_custom_graph()

def make_graph(file):
	f = open(file, 'r')
	lines = f.readlines()
	d = {}
	for line in lines:
		l = line.split(",")
		if float(l[0]) in d.keys():
			d[float(l[0])].append(int(l[1]))
		else:
			d[float(l[0])] = [int(l[1])]
	x = []
	y = []
	keys = d.keys()
	keys.sort()
	for key in keys:
		total = 0
		for val in d[key]:
			total += val
		x.append(key)
		y.append(float(total) / len(d[key]))
	plt.plot(x, y ,'ro-')
	plt.xlabel("Beta")
	plt.ylabel("Iterations")
	plt.show()

def make_custom_graph():
	files = sys.argv[1:]
	file_map = {}
	for file in files:
		f = open(file, 'r')
		lines = f.readlines()
		d = {}
		for line in lines:
			l = line.split(" ")
			if float(l[0]) in d.keys():
				d[float(l[0])].append(int(l[1]))
			else:
				d[float(l[0])] = [int(l[1])]
		x = []
		y = []
		keys = d.keys()
		keys.sort()
		for key in keys:
			total = 0
			for val in d[key]:
				total += val
			x.append(key)
			y.append(float(total) / len(d[key]))
		file_map[file] = (x,y)
	for key in file_map.keys():
		plt.plot(file_map[key][0], file_map[key][1], '-o', label=key.split(':')[0])

	plt.legend(loc=2)
	#plt.axis([0.38, 0.45, -100000, 85000000])
	plt.xlabel("Beta")
	plt.ylabel("Iterations")
	plt.show()


main()