#include "Library.h"
#ifndef DATABASE
#define DATABASE

#include "Database.h"

#endif

using  namespace std;

void lsh_datatype(vector<vector<double>>*, vector<vector<double>>*, vector<int>*, int, int, int, double, vector<int>**, DistanceDatabase<double>*);                 // LSH for Vectors
void lsh_datatype(vector<vector<double*>>*, vector<vector<double*>>*, vector<int>*, int, int, int, double, vector<int>**, DistanceDatabase<double*>*);               // LSH for Curves

template <typename Point> void compute_unassigned(vector<vector<Point>>* ,vector<vector<Point>>* , vector<int>*, int, double**, int**, DistanceDatabase<Point>*);
