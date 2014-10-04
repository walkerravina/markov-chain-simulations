import matplotlib.pyplot as plt
import sys

def main():
	make_graph(sys.argv[1])

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


main()