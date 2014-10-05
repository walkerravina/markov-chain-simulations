import random
import math
import time
import sys

def simulation(n, k):
	"""Run the simulation on a Torus of total size n^2. 
	For each value of beta perform k simulations
	"""
	#vary Beta from small to above critical value
	f = open('results/curie_weiss:' + str(n) + ':' + str(k) + ':time: ' + str(time.time()), 'w')
	beta = 0.01
	while beta <= 2:
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
	#we are on the graph K_n so we represent X and Y simply by lists
	X = [1 for i in range(n)]
	Y = [-1 for i in range(n)]

	#timing information
	iterations = 0
	start = time.clock()
	#+-1 for spins
	global_diff_count = n

	while global_diff_count > 0:
		print(global_diff_count, beta)
		iterations += 1
		#choose a random vertex and spin
		v = random.randint(0,  n - 1)
		new_spin = 1
		if random.random() <= 0.5:
			new_spin = -1

		started_same = (X[v] == Y[v])

		Y_v_spin = Y[v]
		Y_count = 0
		Y_change_count = 0

		for i in range(n):
			if i != v:
				if Y_v_spin != Y[i]:
					Y_count += 1
				if new_spin != Y[i]:
					Y_change_count += 1

		p = min(1, math.exp(-1 * beta * (Y_change_count - Y_count)))

		X_v_spin = X[v]
		X_count = 0
		X_change_count = 0

		for i in range(n):
			if i != v:
				if X_v_spin != X[i]:
					X_count += 1
				if new_spin != X[i]:
					X_change_count += 1

		q = min(1, math.exp(-1 * beta * (X_change_count - X_count)))
		r = random.random()
		
		if r <= p:
			Y[v] = new_spin
		#otherwise leave Y the same
		if r <= q:
			X[v] = new_spin
		#otherwise leave X the same

		#update the count diff
		if started_same and X[v] != Y[v]:
			global_diff_count += 1
		elif not started_same and X[v] == Y[v]:
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