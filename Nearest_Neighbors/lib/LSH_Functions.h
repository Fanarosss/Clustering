#include "Library.h"

using namespace std;

void show_lsh_usage(string name);

template <typename Point> double compute_window(vector<vector<Point>>*);
template <typename Point> void projections(vector<vector<int>>*, vector<vector<Point>>*, vector<double>*, double, int);

void generate_shifts(vector<vector<vector<double>>>*, double, int, int, int);
void compute_hash(vector<int>*, vector<vector<int>>*,int**, int, int, double);
void amplify_hash(vector<int>*, vector<vector<int>>*, int);
