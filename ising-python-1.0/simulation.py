import random
import math
import time
import sys

def simulation(n, k):
	"""Run the simulation on a Torus of total size n^2. 
	For each value of beta perform k simulations
	"""
	#vary Beta from small to above critical value
	f = open('results/ising_model_n:' + str(n) + '_time: ' + str(time.time()), 'w')
	beta = 0.01
	while beta <= .6:
		print("Testing Beta = " + str(beta))
		for i in range(1, k):
			print(i)
			iterations, duration = mix_chains(n, beta)
			f.write(str(beta) + ", " + str(iterations) + ", " + str(duration) + "\n")
		beta += .01
	f.close()

def mix_chains(n, beta):
	"""Run the chains X and Y until they have the same state
	Return the number of moves taken as well as the system time
	"""
	#timing information
	iterations = 0
	start = time.clock()
	#+-1 for spins
	global_diff_count = n * n
	X = [[1 for i in range(n)] for j in range(n)]
	Y = [[-1 for i in range(n)] for j in range(n)]
	shifts = [(1, 0), (-1, 0), (0, 1), (0, -1)] #"shifts" defining the neighbors of a vertex for the torus
	while global_diff_count > 0:
		iterations += 1
		#pick vertex and spin uniformly at random
		v_x, v_y = random.randint(0, n - 1), random.randint(0, n - 1)
		new_spin = 1
		if random.random() <= 0.5:
			new_spin = -1


		started_same = (X[v_x][v_y] == Y[v_x][v_y])

		#there are 4 neighbors, using torus here
		#calc p = Pr(of accepting the change for Y)
		Y_v_spin = Y[v_x][v_y]
		#keep a count for \sigma and \sigma' (\sigma' is where v is assigned the new spin)
		Y_count = 0
		Y_change_count = 0
		for shift in shifts:
			if Y[(v_x + shift[0]) % n][(v_y + shift[1]) % n] != new_spin:
				Y_change_count += 1
			if Y[(v_x + shift[0]) % n][(v_y + shift[1]) % n] != Y_v_spin:
				Y_count += 1
		p = min(1, math.exp(-1 * beta * (Y_change_count - Y_count)))

		#calc q = Pr(of accepting the change for X)
		X_v_spin = X[v_x][v_y]
		X_count = 0
		X_change_count = 0
		for shift in shifts:
			if X[(v_x + shift[0]) % n][(v_y + shift[1]) % n] != new_spin:
				X_change_count += 1
			if X[(v_x + shift[0]) % n][(v_y + shift[1]) % n] != X_v_spin:
				X_count += 1
		q = min(1, math.exp(-1 * beta * (X_change_count - X_count)))
		r = random.random()
		if r <= p:
			Y[v_x][v_y] = new_spin
		#otherwise leave Y the same
		if r <= q:
			X[v_x][v_y] = new_spin
		#otherwise leave X the same
		#
		#update the count diff
		if started_same and X[v_x][v_y] != Y[v_x][v_y]:
			global_diff_count += 1
		elif not started_same and X[v_x][v_y] == Y[v_x][v_y]:
			global_diff_count -= 1

	return (iterations, time.clock() - start)


def main():
	"""
	Main method, read in command line arguments and call simulation
	"""
	if len(sys.argv) <= 2:
		print("Must supply n and k")
	else:
		simulation(int(sys.argv[1]), int(sys.argv[2]))

main()