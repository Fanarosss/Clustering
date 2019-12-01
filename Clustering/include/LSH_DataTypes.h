#include "Library.h"//

using  namespace std;

void lsh_datatype(vector<vector<double>>*, vector<vector<double>>*, int, int, int, vector<int>**);                // LSH for Vectors
void lsh_datatype(vector<vector<double*>>*, vector<vector<double*>>*, int, int, int, vector<int>**);              // LSH for Curves

template <typename Point> void compute_unassigned(vector<vector<Point>>* ,vector<vector<Point>>* , int, double**, int**);
