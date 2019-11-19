#include "Library.h"

using namespace std;

int Read_point_files(vector<vector<int>>*, vector<vector<int>>*, double*, string, string );
int Read_curve_files(vector<vector<double*>>*, vector<vector<double*>>*, string, string);
int Read_curve_files_max_dim(vector<vector<double*>>*, vector<vector<double*>>*, string, string, double);
template <typename Point> double dist(vector<Point>*, vector<Point>*, int, int=1);
double point_dist(double*, double*, int);
double DTW(vector<double*>*, vector<double*>*);
int modulo (int, int);
int moduloMultiplication(int, int, int);
long moduloPower(long, long, long);
void brute_force(vector<vector<int>>*, vector<vector<int>>*, vector<int>*, vector<double>*);
void read_vectors_brute_force_file(string, vector<int>*, vector<double>*);
void curves_brute_force(vector<vector<double*>>*, vector<vector<double*>>*, vector<double>*, vector<double>*, vector<int>*);
void read_curves_brute_force_file(string, vector<double>*, vector<double>*, vector<int>*);
vector<double> arg_min(double**, vector<double>*, double, int);
double min(double, double, double);
