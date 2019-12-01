#include "Initializers.h"
#include "Assigners.h"
#include "Updaters.h"

using namespace std;

template <class Point>
class Cluster {
private:
    Initializer<Point>* initializer;
    Assigner<Point>* assigner;
    Updater<Point>* updater;
    int K;
    int Grids;
    int L;
    int k;
    vector<pair<vector<Point>*, int>> centroids;
    vector<int>** clusters;
public:
    Cluster(int*, string, string, string);
    void fit(vector<vector<Point>>*, DistanceDatabase<Point>*);
    vector<double> silhouette(vector<vector<Point>>*, DistanceDatabase<Point>*);
    double average_distance(int, vector<int>*, DistanceDatabase<Point>*);
    int find_closest_centroid(pair<vector<Point>*, int>, DistanceDatabase<Point>*);
    ~Cluster();
};