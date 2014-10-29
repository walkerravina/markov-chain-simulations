import random
import math
import time
import sys

def simulation(n, k, alpha_low, alpha_high, alpha_step):
	"""Run the simulation on the complete graph of size n. 
	For each value of alpha perform k simulations
	Note that here beta = alpha / n and the critical point should occur at alpha = 1
	"""
	#vary alpha from small to above critical value
	f = open('new_results/curie_weiss:' + str(n) + ':' + str(k) + ':time: ' + str(time.time()), 'w')
	alpha = alpha_low
	while alpha <= alpha_high:
		print("Testing alpha = " + str(alpha))
		for i in range(1, k):
			print(i)
			iterations, duration = mix_chains(n, alpha)
			f.write(str(alpha) + ", " + str(iterations) + ", " + str(duration) + "\n")
		alpha += alpha_step
	f.close()


def mix_chains(n, alpha):
	"""Run the chains X and Y until they have the same state
	Return the number of moves taken as well as the system time
	"""
	#we are on the graph K_n so we represent X and Y by two lists and counts for bookkeeping
	X = [1 for i in range(n)]
	X_pos_total = n
	Y = [-1 for i in range(n)]
	Y_pos_total = 0

	#timing information
	iterations = 0
	start = time.clock()
	#+-1 for spins

	while X_pos_total != Y_pos_total:
		iterations += 1
		#choose a random vertex
		v = random.randint(0,  n - 1)

		#calculate the probability of the vertex being positive with Glauber for both chains
		Y_pos_prob = 0
		Y_pos_count = 0

		if(Y[v] == 1):
			Y_pos_count = Y_pos_total - 1

		spin_sum = Y_pos_count - (n - Y_pos_count - 1)
		Y_pos_prob = math.exp(alpha / n * spin_sum) / (math.exp(alpha / n * spin_sum) + math.exp(-1 * alpha / n * spin_sum))

		X_pos_prob = 0
		X_pos_count = 0

		if(X[v] == 1):
			X_pos_count = X_pos_total - 1

		spin_sum = X_pos_count - (n - X_pos_count - 1)
		X_pos_prob = math.exp(alpha / n * spin_sum) / (math.exp(alpha / n * spin_sum) + math.exp(-1 * alpha / n * spin_sum))

		r = random.random()
		
		if r <= Y_pos_prob:
			if Y[v] == -1:
				Y_pos_total += 1				
			Y[v] = 1
		else:
			if Y[v] == 1:
				Y_pos_total -= 1
			Y[v] = -1
		if r <= X_pos_prob:
			if X[v] == -1:
				X_pos_total += 1
			X[v] = 1
		else:
			if X[v] == 1:
				X_pos_total -= 1
			X[v] = -1

	return (iterations, time.clock() - start)



def main():
	"""
	Main method, read in command line arguments and call simulation
	Command Line Args are:
	size, iterations per alpha, lower alpha, uper alpha, step alpha
	"""
	if len(sys.argv) < 6 :
		print("Must supply n, k, a_low, a_high, a_step space delimited")
	else:
		simulation(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

main()