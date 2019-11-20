#include "Library.h"

using namespace std;

template <typename Point> double compute_window(vector<vector<Point>>*);
template <typename Point> void projections(vector<vector<int>>*, vector<vector<Point>>*, vector<double>*, double, int);

void generate_shifts(vector<vector<vector<double>>>*, double, int, int, int);
void compute_hash(vector<int>*, vector<vector<int>>*,int**, int, int, double);
void amplify_hash(vector<int>*, vector<vector<int>>*, int);

/* Grid Vectorization Functions*/
void Grid_Vectorization(double, int, vector<vector<double*>>*, vector<vector<double*>>*, vector<vector<double>>*, vector<vector<double>>*);
void shift_grid(vector<double>*, int, int);
void hash_curve(vector<vector<double>>*, vector<double*>*, vector<double>*, double, int);
vector<double> arg_min(double**, vector<double>*, double, int);
