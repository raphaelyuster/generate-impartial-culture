// This program generates the F-cultures of size 6 of permutations of {1,2,3,4} and computes the probability of a source vertex (Condorcet winner) in each

#include <iostream>
#include <cmath>


int perm[24][4] = { {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1}, {1,0,2,3}, {1,0,3,2}, {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
					{2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0}, {2,3,0,1}, {2,3,1,0}, {3,0,1,2}, {3,0,2,1}, {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0} };
int maj[4][4];				// maj is the dajacency matrix of the kerel
int generatorIndices[6];	// the indices of the permutations comprising F
int plausible[13][25];		// all plausible values of the coefficients (a,b)  -4 <= a <= 8,  -12 <= b <=12

void printPermutation(int p) { // p is the serial number of the permutation
	int reverse[4];
	for (int i = 0; i < 4; i++)
		reverse[perm[p][i]] = i;
	for (int i = 0; i < 4; i++)
		printf("%d", reverse[i]+1);
	printf(" ");
}

int a, b;
double computeW4probability() { // Compute p(F,W_4) - we shall use maj to host the adjacency matrix of the kernel.
	a = 0, b = 0; // the integer coefficients in the probability expression
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			maj[i][j] = 0;
	for (int t = 0; t < 6; t++)
		for (int i = 0; i < 4; i++)
			for (int j = i + 1; j < 4; j++)
				if (perm[generatorIndices[t]][i] < perm[generatorIndices[t]][j])
					maj[i][j]++;
	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
			if (maj[i][j] > 3)
				maj[i][j] = 1;
			else if (maj[i][j] == 3)
				maj[i][j] = 0;
			else {
				maj[i][j] = 0;
				maj[j][i] = 1;
			}
	int outdeg[4];
	int indeg[4];
	for (int i = 0; i < 4; i++) {
		outdeg[i] = 0;
		indeg[i] = 0;
	}
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			if (maj[i][j] == 1) {
				outdeg[i]++;
				indeg[j]++;
			}
	if (outdeg[0] == 3 || outdeg[1] == 3 || outdeg[2] == 3 || outdeg[3] == 3) {
		a = 4;
		plausible[8][12] = 1; // a=4 b=0
		return 1;
	}
	else if (indeg[0] > 0 && indeg[1] > 0 && indeg[2] > 0 && indeg[3] > 0) {
		plausible[4][12] = 1; // a=0  b=0
		return 0;
	}
	double accumulate = 0;
	for (int u = 0; u < 4; u++) {
		if (indeg[u] > 0)
			continue;
		if (outdeg[u] == 2) {
			accumulate += 0.5;
			a += 2;
		}
		else if (outdeg[u] == 1) {
			// locate v,w
			int v = 0, w = 0;
			for (int i = 0; i < 4; i++)
				if (maj[u][i] == 0 && i != u)
					if (v == 0)
						v = i;
					else
						w = i;
			// compute the correlation between X_{u,v} and X_{u,w}
			int delta = 0;
			for (int t = 0; t < 6; t++)
				if (perm[generatorIndices[t]][u] < perm[generatorIndices[t]][v] && perm[generatorIndices[t]][u] < perm[generatorIndices[t]][w])
					delta++;
			double alpha = 2 * delta / 3.0 - 1;
			accumulate += 0.5 - std::acos(alpha) / (2 * std::acos(-1.0));
			if (delta == 1) {
				a += 2;
				b -= 2;
			}
			else if (delta == 2)
				b += 2;
			else if (delta == 3)
				a += 2;
		}
		else { // outdeg[u] == 0;
			// locate v,w.z
			double val = 0.5;
			a += 2;
			int v = 0, w = 0, z = 0;
			for (int i = 0; i < 4; i++)
				if (maj[u][i] == 0 && i != u)
					if (v == 0)
						v = i;
					else if (w == 0)
						w = i;
					else
						z = i;
			int delta = 0;
			for (int t = 0; t < 6; t++)
				if (perm[generatorIndices[t]][u] < perm[generatorIndices[t]][v] && perm[generatorIndices[t]][u] < perm[generatorIndices[t]][w])
					delta++;
			double alpha = 2 * delta / 3.0 - 1;
			val -= std::acos(alpha) / (4 * std::acos(-1.0));
			if (delta == 0)
				a -= 1;
			else if (delta == 1)
				b -= 1;
			else if (delta == 2) {
				a -= 1;
				b += 1;
			}
			delta = 0;
			for (int t = 0; t < 6; t++)
				if (perm[generatorIndices[t]][u] < perm[generatorIndices[t]][v] && perm[generatorIndices[t]][u] < perm[generatorIndices[t]][z])
					delta++;
			alpha = 2 * delta / 3.0 - 1;
			val -= std::acos(alpha) / (4 * std::acos(-1.0));
			if (delta == 0)
				a -= 1;
			else if (delta == 1)
				b -= 1;
			else if (delta == 2) {
				a -= 1;
				b += 1;
			}
			delta = 0;
			for (int t = 0; t < 6; t++)
				if (perm[generatorIndices[t]][u] < perm[generatorIndices[t]][w] && perm[generatorIndices[t]][u] < perm[generatorIndices[t]][z])
					delta++;
			alpha = 2 * delta / 3.0 - 1;
			val -= std::acos(alpha) / (4 * std::acos(-1.0));
			if (delta == 0)
				a -= 1;
			else if (delta == 1)
				b -= 1;
			else if (delta == 2) {
				a -= 1;
				b += 1;
			}
			accumulate += val;
		}
	}
	plausible[a + 4][b + 12] = 1;
	return accumulate;
}

void main() {
	int count = 0;
	for (int i = 0; i < 13; i++)
		for (int j = 0; j < 25; j++)
			plausible[i][j] = 0;
	for (int i = 0; i < 6; i++)
		generatorIndices[i] = 0;
	for (;;) {
		int i = 1;
		while (generatorIndices[i] == 23)
			generatorIndices[i++] = 0;
		if (i == 6) // exhausted all possibilities
			break;
		generatorIndices[i]++;
		// now compute for F
		double p = computeW4probability();
		if (p < 0.83 && p > 0.815) {
			printf("Permutation set F={ ");
			for(int j=0; j < 6; j++)
				printPermutation(generatorIndices[j]);
			printf("} attains p(F,W_4) = %lf with (a,b) pair = (% d,%d)\n", p,a,b);
			count++;
		}
	}
	printf("Total number of distinct sets attaining p_{W_4} is %d\n", count*24/720);  // multiply by 24 to get all sequences as we always assume that 0123 is the first permutation in the set, and divied by 6! to count sets
	count = 0;
	for (int i = -4; i <= 8; i++)
		for (int j = -12; j <= 12; j++)
			if (plausible[i + 4][j + 12] == 1) {
				printf("A valid (a,b) pair is (%d,%d) giving probability %lf\n", i, j, i / 4.0 + j * std::acos(-1 / 3.0) / (4 * std::acos(-1.0)));
				count++;
			}
	printf("Total number of distinct (a,b) pairs, whence distinct p(F,W_4) probabilities is %d\n", count);
}