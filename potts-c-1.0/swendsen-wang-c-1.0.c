#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000

typedef struct lnode{
	struct lnode* next;
	int val;
} lnode;

void simulation(int n, int k,  double c_low, double c_high, double c_step);
int run_chain(int n, double c, int** spin_assignments, int** visited, lnode*** spin_array, lnode* stk, int* spin_counts);
void llist_add(lnode* head, int val);
int llist_pop(lnode* head);


/**
 * Main method wrapper c
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
 * Run the simulation with the specified parameters. Note that here k=c/n is used in the partition function
 * @param n      The size of the graph (complete graph) n= vertex count
 * @param k      The number of iterations for each c
 * @param c_low  The c to begin with
 * @param c_high The c to end with
 * @param c_step The c to step with 
 */
void simulation(int n, int k, double c_low, double c_high, double c_step){
	int q = 3;
	double c = c_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/swendsen-wang-%d-%d-%d-%f-%f-%f-%d", n, 3, k, c_low, c_high, c_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	int** spin_assignments = malloc(2 * sizeof(int*));
	int** visited = malloc(2 * sizeof(int*));
	lnode*** spin_array = malloc(2 * sizeof(lnode**));
	int* spin_counts = malloc(n * sizeof(int));
	for(int i = 0; i < 2; i++){
		spin_assignments[i] = malloc(n * sizeof(int));
		visited[i] = malloc(n * sizeof(int));
		spin_array[i] = malloc(q * sizeof(lnode*));
		for(int j = 0; j < q; j++){
			spin_array[i][j] = malloc(sizeof(lnode));
			spin_array[i][j]->next = NULL;
			spin_array[i][j]->val = -1; 
			spin_counts[j] = 0;
		}
	}
	lnode* stk = malloc(sizeof(lnode));
	stk->next = NULL;
	stk->val = -1;
	while(c <= c_high){
		for(int i = 0; i < k; i++){
			iterations = run_chain(n, c, spin_assignments, visited, spin_array, stk, spin_counts);
			fprintf(f, "%f %d\n", c, iterations);
			printf("c: %f, k: %d, iterations: %d\n", c, i, iterations);
		}
		c += c_step;
	}
	fclose(f);
}

/**
 * Run the chain until in passes from the type configuration of one dominant spin to all equal distribution of spins
 * @param  n                The number of vertexes
 * @param  c                The value of c in the coupling constant 
 * @param  spin_assignments 2-d array holding the spin assignments for each vertex for the current and next iterations
 * @param  visited          2-d array holding the visited state for each vertex for the current and next iterations	
 * @param  spin_array       2-d array of linked lists holding the vertexes within each spin class
 * @param  stk              stack pointer for the dfs
 * @param  spin_counts      int array for tracking the number of vertexes with each spin on this iteration
 * @return                  The number of iterations required to pass between the two type vectors
 */
int run_chain(int n, double c, int** spin_assignments, int** visited, lnode*** spin_array, lnode* stk, int* spin_counts){
	int q = 3;
	int spin = 0, i = 0, j = 0, k = 0, current = 0, next = 1, equal_spin_count = 0, temp = 0, iterations = 0;
	double p = 1 - exp(-1 * c / n);
	//initialize 
	for(i = 0; i < n; i++){
		if(i < (q - 1) / (double)q * n ){
			spin = 0;
		}
		else if ((q - 1) / (double)q * n <= i && i < (q - 1) / (double)q * n + n / ((q-1) * (double) q)){
			spin = 1;
		}
		else{
			spin = 2;
		}
		spin_assignments[current][i] = spin;
		visited[current][i] = 0;
		visited[next][i] = 0;
		llist_add(spin_array[current][spin], i);
	}

	while(!(equal_spin_count == q)){
		iterations++;
		//for each spin class
		for(i = 0; i < q; i++){
			//create the components for that spin class
			while(spin_array[current][i]->next != NULL){
				j = llist_pop(spin_array[current][i]);
				if(visited[current][j]){
					continue;
				}
				//start the dfs for this new component, pick the spin at random
				spin = arc4random_uniform(q);
				llist_add(stk, j);
				while(stk->next != NULL){
					j = llist_pop(stk);
					if(visited[current][j]){
						continue;
					}
					//mark vertex as visited an update information for the next round
					visited[current][j] = 1;
					visited[next][j] = 0;
					llist_add(spin_array[next][spin], j);
					spin_assignments[next][j] = spin;
					spin_counts[spin]++;
					//for each edge push the neighbor on with the proper probability
					for(k = 0; k < n; k++){
						if(k != j && visited[current][k] == 0 && spin_assignments[current][j] == spin_assignments[current][k] 
							&& ((double)arc4random() / ARC4RANDOM_MAX) <= p){
							llist_add(stk, k);
						}
					}
				}
			}
		}
		equal_spin_count = 0;
		//check if we have passed to type 2 configuration
		for(i = 0; i < q; i++){
			if(spin_counts[i] == n / q){
				equal_spin_count++;
			}
			//reset for next round
			spin_counts[i] = 0;
		}
		//swap current and next
		temp = current;
		current = next;
		next = temp;
	}
	//clear spin array before returning
	for(i = 0; i < q; i++){
		while(spin_array[current][i]->next != NULL){
			llist_pop(spin_array[current][i]);
		}
	}
	return iterations;
}

/**
 * Add the given value to the front of the specified linked list
 * @param head The sentinel head node of the list
 * @param val  The value to add
 */
void llist_add(lnode* head, int val){
	lnode* new_node = malloc(sizeof(lnode));
	new_node->val = val;
	new_node->next = head->next;
	head->next = new_node;
}

/**
 * Pop the head of the given linked list
 * @param  head The sentinel node of the list
 * @return      The value popped 
 */
int llist_pop(lnode* head){
	int ret = head->next->val;
	lnode* ref = head->next;
	head->next = head->next->next;
	free(ref);
	return ret;
}