#include "Library.h"
#include "Database.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
void DistanceDatabase<Point>::calculate_distances(vector<vector<Point>>* cluster_data) {
    /* brute force distances */
    int row, col;
    this->data_size = cluster_data->size();
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
    if (id1 < 0 || id2 < 0) {
        cerr << "Fatal Error: Indexes cannot be less than zero (0)." << endl;
        return -1.0;
    }
    /* convert ids to col and row */
    int row = (id1 > id2) ? id1 : id2;
    int col = (id1 < id2) ? id1 : id2;
    /* out of bounds*/
    if (row >= this->data_size || col >= this->data_size) {
        cerr << "Fatal Error: Index bigger than the Table Size." << endl;
        return -1.0;
    }
    /* same ids = same vectors/curve -> dist = 0 */
    if (row == col) return 0.0;

    return Table[row][col];
}

template <class Point>
DistanceDatabase<Point>::~DistanceDatabase() {
    for (int i = 0; i < this->data_size; i++)
        delete[] this->Table[i];
    delete[] this->Table;
}

template class DistanceDatabase<double>;
template class DistanceDatabase<double*>;