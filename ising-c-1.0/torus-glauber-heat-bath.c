#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000


void simulation(int n, int k, double b_low, double b_high, double b_step);
int mix_chains(int n, double beta);

/**
 * Main method wrapper
 * @param  argc number of args
 * @param  argv arguments array
 * @return      not used
 */
int main(int argc, char *argv[]){
	if (argc != 6){
		printf("Must supply n, k, b_low, b_high, b_step space delimited");
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	double b_low = atof(argv[3]);
	double b_high = atof(argv[4]);
	double b_step = atof(argv[5]);
	simulation(n, k, b_low, b_high, b_step);

}

/**
 * Run the simulation with the specified paramters
 * @param n      The size of the torus
 * @param k      The number of trials to run for each beta
 * @param b_low  The beta to start at
 * @param b_high The beta to end at
 * @param b_step The increment for beta
 */
void simulation(int n, int k, double b_low, double b_high, double b_step){
	double beta = b_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/tours-heat-bath:%d:%d:%f:%f:%f:%d", n, k, b_low, b_high, b_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	while(beta <= b_high){
		for(int i = 0; i < k; i++){
			iterations = mix_chains(n, beta);
			fprintf(f, "%f %d\n", beta, iterations);
			printf("beta: %f, k: %d, iterations: %d\n", beta, i, iterations);
		}
		beta += b_step;
	}
	fclose(f);
}

/**
 * Run the chains X and Y until they couple
 * @param  n    The size of the 2D torus
 * @param  beta The value for beta for the partition function
 * @param  X    Pointer for the X chain
 * @param  Y    Pointer for the Y chain
 * @return      The iterations required for coupling
 */
int mix_chains(int n, double beta){
	int X[n][n];
	int Y[n][n];
	unsigned long long iterations = 0;
	/**
	 * use +1 / -1 for spins, stat X at all + and Y at all -
	 */
	
	int global_diff_count = n * n;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			X[i][j] = 1;
			Y[i][j] = -1;
		}
	}

	int v_x, v_y, started_same, local_spin_sum,
		v_x_new, v_y_new;
	double Y_pos_prob, X_pos_prob, r;

	while(global_diff_count > 0){
		iterations += 1;
		v_x = arc4random_uniform(n);
		v_y = arc4random_uniform(n);

		started_same = X[v_x][v_y] - Y[v_x][v_y];

		local_spin_sum = 0;
		//sum the spins of the neighbors observing wrap around in the Y chain
		v_x_new = (v_x + 1) % n;
		v_y_new = v_y;
		local_spin_sum += Y[v_x_new][v_y_new];
		v_x_new = (v_x - 1) % n;
		if(v_x_new < 0){
			v_x_new += n;
		}
		local_spin_sum += Y[v_x_new][v_y_new];

		v_x_new = v_x;
		v_y_new = (v_y + 1) % n;
		local_spin_sum += Y[v_x_new][v_y_new];
		v_y_new = (v_y - 1) % n;
		if(v_y_new < 0){
			v_y_new += n;
		}
		local_spin_sum += Y[v_x_new][v_y_new];

		//calculate the probability using glauber dynamics for Y chain
		Y_pos_prob = exp(beta * local_spin_sum) / (exp(beta * local_spin_sum ) + exp(-1 * beta * local_spin_sum) );


		local_spin_sum = 0;

		//sum the spins of the neighbors observing wrap around in the X chain
		v_x_new = (v_x + 1) % n;
		v_y_new = v_y;
		local_spin_sum += X[v_x_new][v_y_new];
		v_x_new = (v_x - 1) % n;
		if(v_x_new < 0){
			v_x_new += n;
		}
		local_spin_sum += X[v_x_new][v_y_new];

		v_x_new = v_x;
		v_y_new = (v_y + 1) % n;
		local_spin_sum += X[v_x_new][v_y_new];
		v_y_new = (v_y - 1) % n;
		if(v_y_new < 0){
			v_y_new += n;
		}
		local_spin_sum += X[v_x_new][v_y_new];


		//calculate the probability using glauber dynamics for Y chain
		X_pos_prob = exp(beta * local_spin_sum) / (exp(beta * local_spin_sum ) + exp(-1 * beta * local_spin_sum));


		r = ((double)arc4random() / ARC4RANDOM_MAX);
		
		//update the chains with the proper probabilities 
		if (r <= Y_pos_prob){
			Y[v_x][v_y] = 1;
		}
		else{
			Y[v_x][v_y] = -1;
		}
		if (r <= X_pos_prob){
			X[v_x][v_y] = 1;
		}
		else{
			X[v_x][v_y] = -1;
		}

		//keep track of how many vertexes are different
		if (started_same == 0 && X[v_x][v_y] != Y[v_x][v_y]){
			global_diff_count += 1;
		}
		else if(started_same != 0 && X[v_x][v_y] == Y[v_x][v_y]){
			global_diff_count -= 1;
		}		
	}


	return iterations;

}