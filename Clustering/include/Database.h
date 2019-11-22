#include "Library.h"

using namespace std;

template <class Point>
class DistanceDatabase {
private:
    double** Table;
public:
    void calculate_distances(vector<vector<Point>>*);
    double get_distance(int, int);
};
