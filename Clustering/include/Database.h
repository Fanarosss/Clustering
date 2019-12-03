#include "Library.h"

using namespace std;

template <class Point>
class DistanceDatabase {
private:
    double** Table;
    int data_size;
public:
    void calculate_distances(vector<vector<Point>>*);
    double get_distance(int, int);
    ~DistanceDatabase();
};
