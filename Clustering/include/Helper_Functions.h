#include "Library.h"

using namespace std;

template <typename Point> double min_distance(int, vector<int>*, vector<vector<Point>>*);

void show_cluster_usage(string name);
int Read_input_file(string input);
int Read_files(vector<vector<double>>*, int*, string, string, vector<string>*);
int Read_files(vector<vector<double*>>*, int*, string, string, vector<string>*);
double point_dist(double*, double*, int=2);
double DTW(vector<double*>*, vector<double*>*);
void DTW_pairs(vector<double*>*, vector<double*>*, vector<pair<int,int>>*);
double dist(vector<int>* P1, vector<int>* P2, int, int=1);
double dist(vector<double>* P1, vector<double>* P2, int, int=1);
double dist(vector<double*>* P1, vector<double*>* P2, int, int=1);
double min(double, double, double);
void normalize(vector<double>*);
double Sum(int, int, vector<double>*, int);
int modulo (int, int);
int moduloMultiplication(int, int, int);
long moduloPower(long, long, long);
