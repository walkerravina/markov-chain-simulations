#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000


void simulation(int n, int k, double a_low, double a_high, double a_step);
int mix_chains(int n, double beta);

/**
 * Main method wrapper
 * @param  argc the number of arguments
 * @param  argv the arguments array
 * @return      not used
 */
int main(int argc, char *argv[]){
	if (argc != 6){
		printf("Must supply n, k, a_low, a_high, a_step space delimited");
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	double a_low = atof(argv[3]);
	double a_high = atof(argv[4]);
	double a_step = atof(argv[5]);
	simulation(n, k, a_low, a_high, a_step);

}

/**
 * Run the simulation with the specified parameters. Note that here alpha/n is used in the partition function,
 * this is also referred to as beta/n
 * @param n      The size of the graph (complete graph) n= vertex count
 * @param k      The number of iterations for each alpha
 * @param a_low  The alpha to begin with
 * @param a_high The alpha to end with
 * @param a_step The alpha to step with 
 */
void simulation(int n, int k, double a_low, double a_high, double a_step){
	double alpha = a_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/curie-weiss-heat-bath:%d:%d:%f:%f:%f:%d", n, k, a_low, a_high, a_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	while(alpha <= a_high){
		for(int i = 0; i < k; i++){
			iterations = mix_chains(n, alpha);
			fprintf(f, "%f %d\n", alpha, iterations);
			printf("alpha: %f, k: %d, iterations: %d\n", alpha, i, iterations);
		}
		alpha += a_step;
	}
	fclose(f);
}

/**
 * Run an all positive and all negative starting chains until they couple
 * and report the required number of iterations
 * @param  n     The size of the chains
 * @param  alpha The alpha for the partition function
 * @return       The iterations needed for mixing
 */
int mix_chains(int n, double alpha){
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

	while (X_pos_total != Y_pos_total){

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
		
		Y_pos_prob = exp(alpha / n * spin_sum) / (exp(alpha / n * spin_sum) + exp(-1 * alpha / n * spin_sum));


		if(X[v] == 1){
			spin_sum = X_pos_total - 1 - (n - X_pos_total);		
		}
		else{
			spin_sum = X_pos_total - (n - X_pos_total - 1);
		}

		X_pos_prob = exp(alpha / n * spin_sum) / (exp(alpha / n * spin_sum) + exp(-1 * alpha / n * spin_sum));

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