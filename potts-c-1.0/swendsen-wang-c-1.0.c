#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define ARC4RANDOM_MAX      0x100000000

typedef struct lnode{
	struct lnode* next;
	int val;
} lnode;

void simulation(int n, int k, int q, double c_low, double c_high, double c_step);
int run_chain(int n, int q, double c, int* spin_assignments, lnode** spin_array, int* spin_counts, lnode** edge_set);
void llist_add(lnode* head, int val);
void edge_set_add(lnode** edge_set, int i, int j);
int llist_pop(lnode* head);
int edge_set_find(lnode** edge_set, int i, int j);
void clean_data_structures(lnode** edge_set, lnode** spin_array, int* spin_counts, int n, int q);


/**
 * Main method wrapper c
 * @param  argc the number of arguments
 * @param  argv the arguments array
 * @return      not used
 */
int main(int argc, char *argv[]){
	if (argc != 7){
		printf("Must supply n, k, q, c_low, c_high, c_step space delimited");
	}
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	int q = atoi(argv[3]);
	double c_low = atof(argv[4]);
	double c_high = atof(argv[5]);
	double c_step = atof(argv[6]);
	simulation(n, k, q, c_low, c_high, c_step);

}

/**
 * Run the simulation with the specified parameters. Note that here k=c/n is used in the partition function
 * @param n      The size of the graph (complete graph) n= vertex count
 * @param q  	 The number of spins
 * @param k      The number of iterations for each c
 * @param c_low  The c to begin with
 * @param c_high The c to end with
 * @param c_step The c to step with 
 */
void simulation(int n, int k, int q, double c_low, double c_high, double c_step){
	double c = c_low;
	int iterations;
	char file_name[100];
	sprintf(file_name, "results/swendsen-wang-%d-%d-%d-%f-%f-%f-%d", n, q, k, c_low, c_high, c_step, (unsigned)time(NULL));
	FILE *f = fopen(file_name, "w");
	if(NULL == f){
		printf("Error opening results file");
		exit(1);
	}
	//initialize an array to keep track of the number of vertexes with each spin
	//and an array to keep track of the count of vertexes with each spin
	lnode** spin_array = malloc(q * sizeof(lnode));
	int* spin_counts = malloc(q * sizeof(int));
	for(int i = 0; i < q; i++){
		//use a sentinel for llists
		spin_array[i] = malloc(sizeof(lnode));
		spin_array[i]->next = NULL;
		spin_array[i]->val = -1;
		spin_counts[i] = 0;
	}
	//edge set serves as a hash set for tracking which edges have been removed in phase one of swendsen wang
	//there are n buckets, each vertex i has a bucket with a list of vertexes such that edge (i, j) has been removed
	//for each vertex j in the bucket. When an edge (a, b) is added to the set the bucket which is used is for the lower id vertex
	//edges should be checked for membership in a similar way
	lnode** edge_set = malloc(n * sizeof(lnode));
	int* spin_assignments = malloc(n * sizeof(int));
	for(int i = 0; i < n; i++){
		edge_set[i] = malloc(sizeof(lnode));
		//use a sentinel
		edge_set[i]->val = 0;
		edge_set[i]->next = NULL;
		spin_assignments[i] = 0;
	}
	while(c <= c_high){
		for(int i = 0; i < k; i++){
			clean_data_structures(edge_set, spin_array, spin_counts, n, q);
			iterations = run_chain(n, q, c, spin_assignments, spin_array, spin_counts, edge_set);
			fprintf(f, "%f %d\n", c, iterations);
			printf("c: %f, k: %d, iterations: %d\n", c, i, iterations);
		}
		c += c_step;
	}
	fclose(f);
}

/**
 * Run a chain starting at the type vector T such that T[i] = 1/q for all i
 * that is the equal distribution type vector. We run the swendsen wang process
 * until the type vector becomes T' such that T'[k] = (q-1)/2 and T' [i!=k] = 1/(q(q-1))
 * that is the type vector where one spin dominates
 * @param  n    The number of vertexes defining the complete graph
 * @param  q    The number of spins
 * @param  c The c in the k=c/n value for the coupling constant
 * @param  spin_array Array of linked lists for keeping track of vertexes with each spin
 * @param  spin_counts Array of integers for keeping track of count of vertexes with each spin
 * @param  edge_set Hash table for keeping track of removed edges
 * @return      The number of iterations required to reach the desired Type vector
 */
int run_chain(int n, int q, double c, int* spin_assignments, lnode** spin_array, int* spin_counts, lnode** edge_set){
	int spin = 0; //spins range from 0 to q-1 for easy indexing
	int i = 0, j = 0, k = 0;
	//initialize the system to the first type vector
	for(i = 0; i < n; i++){
		spin = i % 3;
		spin_counts[spin]++;
		spin_assignments[i] = spin;
		llist_add(spin_array[spin], i);
	}
	int dom_spin_count = 0, other_spin_count = 0, iterations = 0;
	double p = exp(-1 * c / n);//probability of removing an edge
	lnode* v1 = NULL; lnode* curr = NULL; lnode* stk = malloc(sizeof(lnode));
	//stk is for dfs, again we use a sentinel node
	stk->val = -1;
	stk->next = NULL;
	while(!(dom_spin_count == 1 && other_spin_count == q - 1)){
		iterations++;
		//remove edges
		//for each spin group
		for(i = 0; i < q; i++){
			//for each unordered pair of vertexes (i.e each edge) in that group
			//remove that edge with the proper probability, we shrink the list as we go
			while(spin_array[i]->next != NULL){
				v1 = spin_array[i]->next;
				curr = v1->next;
				while(curr != NULL){
					if(((double)arc4random() / ARC4RANDOM_MAX) <= p){
						edge_set_add(edge_set, v1->val, curr->val);
					}
					curr = curr->next;
				}
				spin_array[i]->next = spin_array[i]->next->next;
				free(v1);
			}
			//reset spin counts
			spin_counts[i] = 0;
		}
		//perform DFS to randomly assign a spin to each connected component
		//edge set doubles as a visited list with the sentinel nodes values tracking
		//whether the node has been visited
		for(i = 0; i < n; i++){
			if(edge_set[i]->val){
				continue;
			}
			//chose a random spin for this connected component
			spin = arc4random_uniform(q);
			llist_add(stk, i);
			while(stk->next != NULL){
				j = llist_pop(stk);
				if(edge_set[j]->val){
					continue;
				}
				edge_set[j]->val = 1;
				llist_add(spin_array[spin], j);
				spin_counts[spin]++;
				for(k = 0; k < n; k++){
					//if this node has already been visited, or the edge is not actually there then skip
					//otherwise push onto the stk
					if(k == j ||  edge_set[k]->val == 1 || spin_assignments[k] != spin_assignments[j] || edge_set_find(edge_set, j, k) == 1 ){
						continue;
					}
					spin_assignments[j] = spin;
					llist_add(stk, k);
				}
			}
		}
		//reset visited list for next iteration
		for(i = 0; i < n; i++){
			if(iterations == 1){
				curr = edge_set[i];
				printf("slot %d: \n",i );
				while(curr != NULL){
					printf("%d, \n", curr->val);
					curr = curr->next;
				}
				printf("\n");
			}
			edge_set[i]->val = 0;
		}
		dom_spin_count = 0; other_spin_count = 0;
		//check if we have passed to type 2 configuration
		for(i = 0; i < q; i++){
			if(spin_counts[i] == n * (q-1) / q){
				dom_spin_count++;
			}
			else if(spin_counts[i] == n / (q*(q-1)) ){
				other_spin_count++;
			}
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

/**
 * Add the given (i, j) edge to our edge set
 * Edge is bucketed under a = min(i,j) with b = max(i,j) added to the linked list in that bucket
 * @param edge_set Pointer to the edge set
 * @param i        vertex 1
 * @param j        vertex 2
 */
void edge_set_add(lnode** edge_set, int i, int j){
	lnode* head;
	lnode* new_node = malloc(sizeof(lnode));
	int lower = j, higher = i;
	if(i < j){
		lower = i;
		higher = j;
	}
	head = edge_set[lower];
	new_node->val = higher;
	new_node->next = head->next;
	head->next = new_node;
}

/**
 * Determine whether the specified edge is contained in our edge set
 * after being checked an edge is removed since it will only need to be checked once
 * @param  edge_set Pointer to the edge set
 * @param  i        vertex 1
 * @param  j        vertex 2
 * @return          1 if edge is in the set, 0 otherwise
 */
int edge_set_find(lnode** edge_set, int i, int j){
	int lower = j, higher = i;
	if(i < j){
		lower = i;
		higher = j;
	}
	lnode* prev = edge_set[lower];
	lnode* curr = edge_set[lower]->next;
	while(curr != NULL){
		//1st time reading this edge
		if(curr->val == higher){
			prev->next = curr->next;
			free(curr);
			return 1;
		}
		prev = curr;
		curr = curr->next;
	}
	return 0;
}

/**
 * Clean all of the data structures for the next run
 * @param edge_set    the edge_set
 * @param spin_array  the spin array
 * @param spin_counts the spin counts array
 * @param n           the number of vertexes
 * @param q           the number of spins
 */
void clean_data_structures(lnode** edge_set, lnode** spin_array, int* spin_counts, int n, int q){
	for(int i = 0; i < q; i++){
		spin_array[i]->next = NULL;
		spin_array[i]->val = -1;
		spin_counts[i] = 0;
	}
	for(int i = 0; i < n; i++){
		edge_set[i]->val = 0;
		edge_set[i]->next = NULL;
	}
}