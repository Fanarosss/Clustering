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
    vector<int>* centroids;
    vector<int>** clusters;
public:
    Cluster(int*, string, string, string);
    void fit(vector<vector<Point>>*, DistanceDatabase<Point>*);
    vector<double> silhouette(vector<vector<Point>>*);
    double average_distance(int, vector<int>*, vector<vector<Point>>*);
    int find_closest_centroid(int, vector<int>*, vector<vector<Point>>*);
    ~Cluster();
};