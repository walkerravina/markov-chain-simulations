import random
import math
import time
import sys

def simulation(n, k, b_low, b_high, b_step):
	"""Run the simulation on a Torus of total size n^2. 
	For each value of beta perform k simulations
	"""
	#vary Beta from small to above critical value
	file_name = 'results/torus-heat-bath:' + str(n) + ':k:' + str(k) + ':beta_low:' + str(b_low) + ':beta_high:'
	file_name += str(b_high) + ':beta_step:' + str(b_step) + ':time:' + str(time.time())
	f = open(file_name, 'w')
	beta = b_low
	while beta <= b_high:
		print("Testing Beta = " + str(beta))
		for i in range(0, k):
			print(i)
			iterations, duration = mix_chains(n, beta)
			f.write(str(beta) + ", " + str(iterations) + ", " + str(duration) + "\n")
		beta += b_step
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
		#pick vertex uniformly at random
		v_x, v_y = random.randint(0, n - 1), random.randint(0, n - 1)

		started_same = (X[v_x][v_y] == Y[v_x][v_y])

		#there are 4 neighbors, using torus here
		#determine P(change to +) Y
		Y_pos_prob = 0
		local_spin_sum = 0
		for shift in shifts:
			local_spin_sum += Y[(v_x + shift[0]) % n][(v_y + shift[1]) % n]

		Y_pos_prob = math.exp(beta * local_spin_sum) / (math.exp(beta * local_spin_sum ) + math.exp(-1 * beta * local_spin_sum) )

		#determine P(change to +) for X
		X_pos_prob = 0
		local_spin_sum = 0

		for shift in shifts:
			local_spin_sum += X[(v_x + shift[0]) % n][(v_y + shift[1]) % n]

		X_pos_prob = math.exp(beta * local_spin_sum) / (math.exp(beta * local_spin_sum ) + math.exp(-1 * beta * local_spin_sum) )
		r = random.random()
		
		if r <= Y_pos_prob:
			Y[v_x][v_y] = 1
		else:
			Y[v_x][v_y] = -1
		if r <= X_pos_prob:
			X[v_x][v_y] = 1
		else:
			X[v_x][v_y] = -1

		#update the count diff
		if started_same and X[v_x][v_y] != Y[v_x][v_y]:
			global_diff_count += 1
		elif not started_same and X[v_x][v_y] == Y[v_x][v_y]:
			global_diff_count -= 1

	return (iterations, time.clock() - start)


def main():
	"""
	Main method, read in command line arguments and call simulation
	Command Line Args are:
	size, iterations per beta, lower beta, uper beta, step beta
	"""
	if len(sys.argv) < 6 :
		print("Must supply n, k, b_low, b_high, b_step space delimited")
	else:
		simulation(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

main()