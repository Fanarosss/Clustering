#include "Library.h"

using namespace std;

void show_grid_lsh_usage(string name);
void show_grid_bhc_usage(string name);
void shift_grid(vector<double>*, int, int);
void hash_curve(vector<vector<double>>*, vector<double*>*, vector<double>*, double, int);
void find_diagonal(vector<vector<int>>*, int, int);
void find_traversals(vector<vector<int>>*, int, int, int, int, vector<vector<int>>*, vector<vector<vector<int>>>*);
void Relevant_Traversals(vector<vector<vector<int>>>*, int, int);
