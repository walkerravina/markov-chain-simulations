#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000


void simulation(int n, int k, double lambda_low, double lambda_high, double lambda_step);
int mix_chains(int n, double beta);

/**
 * Main method wrapper
 * @param  argc the number of arguments
 * @param  argv the arguments array
 * @return      not used
 */
int main(int argc, char *argv[]){
	if (argc != 6){
		printf("Must supply n, k, lambda_low, lambda_high, lambda_step space delimited");
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	double lambda_low = atof(argv[3]);
	double lambda_high = atof(argv[4]);
	double lambda_step = atof(argv[5]);
	simulation(n, k, lambda_low, lambda_high, lambda_step);

}

/**
 * Run the simulation with the specified parameters.
 * @param n      The size of the graph (complete graph) n= vertex count
 * @param k      The number of iterations for each lambda
 * @param a_low  The lambda to begin with
 * @param a_high The lambda to end with
 * @param a_step The lambda to step with 
 */
void simulation(int n, int k, double lambda_low, double lambda_high, double lambda_step){
	double lambda = lambda_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/independent-set-heat-bath:%d:%d:%f:%f:%f:%d", n, k, lambda_low, lambda_high, lambda_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	while(lambda <= lambda_high){
		for(int i = 0; i < k; i++){
			iterations = mix_chains(n, lambda);
			fprintf(f, "%f %d\n", lambda, iterations);
			printf("lambda: %f, k: %d, iterations: %d\n", lambda, i, iterations);
		}
		lambda += lambda_step;
	}
	fclose(f);
}

/**
 * Run even and odd occupied starting chains until they couple
 * and report the required number of iterations
 * @param  n     The size of the chains
 * @param  lambda The lambda for the partition function
 * @return       The iterations needed for mixing
 */
int mix_chains(int n, double lambda){
	//we are on the graph K_n so we represent X and Y by two lists and counts for bookkeeping
	int X[n][n];
	int Y[n][n];

	int global_diff_count = n * n;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if((i +j) % 2 == 0){
				X[i][j] = 1;
				Y[i][j] = 0;
			}
			else{
				X[i][j] = 0;
				Y[i][j] = 1;
			}
		}
	}


	unsigned long long iterations = 0;
	// +1 / -1 for spins

	int v_x, v_y, flag, v_x_new, v_y_new, started_same;
	double occupation_prob, r;

	while (global_diff_count > 0){

		iterations += 1;
		v_x = arc4random_uniform(n);
		v_y = arc4random_uniform(n);
		started_same = X[v_x][v_y] - Y[v_x][v_y];

		occupation_prob = lambda / (lambda + 1);

		r = ((double)arc4random() / ARC4RANDOM_MAX);
		
		if (r <= occupation_prob){
			//check the neighbors to ensure valid IS
			flag = 1;
			v_x_new = (v_x + 1) % n;
			if(Y[v_x_new][v_y] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_y_new = (v_y + 1) % n;
			if(Y[v_x][v_y_new] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_x_new = v_x - 1;
			if(v_x_new < 0){
				v_x_new += n;
			}
			if(Y[v_x_new][v_y] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_y_new = v_y - 1;
			if(v_y_new < 0){
				v_y_new += n;
			}
			if(Y[v_x][v_y_new] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			if(flag){
				Y[v_x][v_y] = 1;
			}
		}
		else{
			Y[v_x][v_y] = 0;	
		}
		if (r <= occupation_prob){
			//check the neighbors to ensure valid IS
			flag = 1;
			v_x_new = (v_x + 1) % n;
			if(X[v_x_new][v_y] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_y_new = (v_y + 1) % n;
			if(X[v_x][v_y_new] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_x_new = v_x - 1;
			if(v_y_new < 0){
				v_x_new += n;
			}
			if(X[v_x_new][v_y] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			v_y_new = v_y - 1;
			if(v_y_new < 0){
				v_y_new += n;
			}
			if(X[v_x][v_y_new] != 1){
				flag &= 1;
			}
			else{
				flag &= 0;
			}
			if(flag){
				X[v_x][v_y] = 1;
			}		
		}
		else{
			X[v_x][v_y] = 0;		
		}
		if(started_same == 0 && X[v_x][v_y] != Y[v_x][v_y]){
			global_diff_count += 1;
		}
		else if(started_same != 0 && X[v_x][v_y] == Y[v_x][v_y]){
			global_diff_count -= 1;
		}
	}
	return iterations;
}