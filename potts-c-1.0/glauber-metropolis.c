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
 * as tge coupling constant
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
	sprintf(file_name, "results/curie-weiss-heat-bath:%d:%d:%f:%f:%f:%d", n, k, c_low, c_high, c_step, (unsigned)time(NULL));
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
 * Run q chains each starting at a configuration in which all of the vertexes are the same color
 * @param  n     The size of the chains
 * @param  c The c for the partition function
 * @return       The iterations needed for mixing
 */
int mix_chains(int n, double c){
	//we are on the graph K_n so we represent X and Y by two lists and counts for bookkeeping
	int X[n];
	int Y[n];

	for(int i = 0; i < n; i++){
		X[i] = 1;
		Y[i] = -1;
	}

	int X_pos_total = n;
	int Y_pos_total = 0;
	int global_diff_count = n;

	unsigned long long iterations = 0;
	// +1 / -1 for spins

	int v, spin_sum, started_same;
	double Y_pos_prob, X_pos_prob, r;

	while (global_diff_count > 0){

		iterations += 1;
		v = arc4random_uniform(n);
		started_same = X[v] - Y[v];

		//calculate the probability of the vertex being positive with Glauber for both chains

		if(Y[v] == 1){
			spin_sum = Y_pos_total - 1 - (n - Y_pos_total);
		}
		else{
			spin_sum = Y_pos_total - (n - Y_pos_total - 1);
		}
		
		Y_pos_prob = exp(c / n * spin_sum) / (exp(c / n * spin_sum) + exp(-1 * c / n * spin_sum));


		if(X[v] == 1){
			spin_sum = X_pos_total - 1 - (n - X_pos_total);		
		}
		else{
			spin_sum = X_pos_total - (n - X_pos_total - 1);
		}

		X_pos_prob = exp(c / n * spin_sum) / (exp(c / n * spin_sum) + exp(-1 * c / n * spin_sum));

		r = ((double)arc4random() / ARC4RANDOM_MAX);
		
		if (r <= Y_pos_prob){
			if(Y[v] != 1){
				Y_pos_total += 1;
			}	
			Y[v] = 1;
		}
		else{
			if(Y[v] == 1){
				Y_pos_total -= 1;
			}
			Y[v] = -1;		
		}

		if (r <= X_pos_prob){
			if(X[v] != 1){
				X_pos_total += 1;
			}
			X[v] = 1;			
		}
		else{
			if (X[v] == 1){
				X_pos_total -= 1;
			}
			X[v] = -1;			
		}
		if(started_same == 0 && X[v] != Y[v]){
			global_diff_count += 1;
		}
		else if(started_same != 0 && X[v] == Y[v]){
			global_diff_count -= 1;
		}
	}

	return iterations;
}