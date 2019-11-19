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
    vector<int>* centroids;
    vector<int>** clusters;
public:
    Cluster(int, string, string, string);
    void fit(vector<vector<Point>>*);
    ~Cluster();
};
