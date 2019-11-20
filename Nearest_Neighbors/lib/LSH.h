#include "HashTable.h"

using namespace std;

template <class Point>
class LSH {
private:
    /* dataset */
    vector<vector<Point>>* dataset;
    /* Vector containing shifts of size(l,k,d) */
    vector<vector<vector<double>>> s;
    /* Amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* Vector containing projections of data */
    vector<vector<int>> a_projects;
    /* L Hash Tables */
    HashTable <Point> **MyHashTable;
    /* Size of hash Table */
    int TableSize;
    /* array for hash computations */
    int * power;
    /* parameters LSH */
    int k;
    int L;
    Point w;
    int d;
    /* end of params */
public:
    LSH(int, int, Point);
    void fit(vector<vector<Point>>*);
    void evaluate(vector<vector<Point>>*, double, vector<vector<int>>*, Point**, double**, int**);
    ~LSH();
};

void Grid_Vectorization(double, int, vector<vector<double*>>*, vector<vector<double*>>*, vector<vector<double>>*, vector<vector<double>>*);