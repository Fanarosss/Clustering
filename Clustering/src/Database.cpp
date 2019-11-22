#include "Library.h"
#include "Database.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
void DistanceDatabase<Point>::calculate_distances(vector<vector<Point>>* cluster_data) {
    /* brute force distances */
    int row, col;
    int data_size = cluster_data->size();
    int dimensions = (*cluster_data)[0].size();
    this->Table = new double* [data_size];
    /* init array */
    for (int i = 0; i < data_size; i++) {
        Table[i] = new double [data_size];
        for (int j = 0; j < data_size; j++) {
            Table[i][j] = -1;
        }
    }
    /* calculate dists */
    for (int i = 0; i < data_size; i++) {
        for (int j = 0; j < data_size; j++) {
            row = (i > j) ? i : j;
            col = (i < j) ? i : j;
            if (row != col) {
                if (Table[row][col] == -1)
                    Table[row][col] = dist(&(*cluster_data)[row], &(*cluster_data)[col], dimensions);
            } else {
                Table[row][col] = 0;
            }
        }
    }
    /* end of brute force */
}

template <class Point>
double DistanceDatabase<Point>::get_distance(int id1, int id2) {
    int row = (id1 > id2) ? id1 : id2;
    int col = (id1 < id2) ? id1 : id2;
    return Table[row][col];
}

template class DistanceDatabase<int>;
template class DistanceDatabase<double*>;