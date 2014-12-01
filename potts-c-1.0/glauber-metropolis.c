#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000


void simulation(int n, int k, double c_low, double c_high, double c_step);
int mix_chains(int n, double c);

/**
 * Main method wrapper
 * @param  argc the number of arguments
 * @param  argv the arguments array
 * @return      not used
 */
int main(int argc, char *argv[]){
	if (argc != 6){
		printf("Must supply n, k, c_low, c_high, c_step space delimited");
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	double c_low = atof(argv[3]);
	double c_high = atof(argv[4]);
	double c_step = atof(argv[5]);
	simulation(n, k, c_low, c_high, c_step);

}

/**
 * Run the simulation with the specified parameters. Note that here c/n is used in the partition function,
 * as the coupling constant
 * @param n      The size of the graph (complete graph) n= vertex count
 * @param k      The number of iterations for each c
 * @param c_low  The c to begin with
 * @param c_high The c to end with
 * @param c_step The c to step with 
 */
void simulation(int n, int k, double c_low, double c_high, double c_step){
	double c = c_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/glauber-metropolis-%d-%d-%f-%f-%f-%d", n, k, c_low, c_high, c_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	while(c <= c_high){
		for(int i = 0; i < k; i++){
			iterations = mix_chains(n, c);
			fprintf(f, "%f %d\n", c, iterations);
			printf("c: %f, k: %d, iterations: %d\n", c, i, iterations);
		}
		c += c_step;
	}
	fclose(f);
}

/**
 * Run 3 chains (q=3) each starting at a configuration in which all of the vertexes have the same spin
 * @param  n     The size of the chains
 * @param  c The c for the partition function
 * @return       The iterations needed for mixing
 */
int mix_chains(int n, double c){
	//we are on the graph K_n so we represent X, Y, Z by lists and type vectors
	int X[n];
	int Y[n];
	int Z[n];
	int X_type[3];
	int Y_type[3];
	int Z_type[3];

	for(int i = 0; i < n; i++){
		X[i] = 0;
		Y[i] = 1;
		Z[i] = 2;
	}
	for(int i = 0; i < 3; i++){
		X_type[i] = 0;
		Y_type[i] = 0;
		Z_type[i] = 0;
	}
	X_type[0] = n;
	Y_type[1] = n;
	Z_type[2] = n;

	unsigned long long iterations = 0;
	// +1 / -1 for spins

	double X_prob, Y_prob, Z_prob, r;
	int not_done = 1, new_spin = 0, old_spin = 0, v= 0;
	//run the chains until the type vectors match up
	while (not_done){
		not_done = 0;
		for(int i = 0; i < 3; i++){
			if(X[i] == Y[i] && X[i] == Z[i]){
				continue;
			}
			else{
				not_done = not_done | 1;
			}
		}
		iterations += 1;
		v = arc4random_uniform(n);
		new_spin = arc4random_uniform(3);

		//make the move with the proper probability in each chain according to the metropolis rule		

		r = ((double)arc4random() / ARC4RANDOM_MAX);

		old_spin = X[v];
		X_prob = X_type[new_spin] - (X_type[old_spin] - 1);
		X_prob = exp(c / n * X_prob);
		if(r <= X_prob){
			X[v] = new_spin;
			X_type[old_spin]--;
			X_type[new_spin]++;
		}
				
		old_spin = Y[v];
		Y_prob = Y_type[new_spin] - (Y_type[old_spin] - 1);
		Y_prob = exp(c / n * Y_prob);
		if(r <= Y_prob){
			Y[v] = new_spin;
			Y_type[old_spin]--;
			Y_type[new_spin]++;
		}	

		old_spin = Z[v];
		Z_prob = Z_type[new_spin] - (Z_type[old_spin] - 1);
		Z_prob = exp(c / n * Z_prob);
		if(r <= Z_prob){
			Z[v] = new_spin;
			Z_type[old_spin]--;
			Z_type[new_spin]++;
		}
	}

	return iterations;
}