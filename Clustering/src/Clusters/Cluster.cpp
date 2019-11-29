#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
Cluster <Point>::Cluster(int* cluster_conf, string Initializer, string Assigner, string Updater) {
    /* number of clusters*/
    this->K = cluster_conf[0];
    /* number of grids*/
    this->Grids = cluster_conf[1];
    /* number of vector hash tables*/
    this->L = cluster_conf[2];
    /* number of vector hash functions*/
    this->k = cluster_conf[3];

    /* Initializer */
    cout << "------Configuration------" << endl;
    if (Initializer == "Random Selection") {
        this->initializer = new Random_Selection<Point>(this->K);
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(this->K);
    } else {
        cerr << "Unknown Initializer";
    }
    cout << '\t' << "Initializer: " << initializer->get_name() << endl;
    /* Assigner */
    if (Assigner == "Lloyd's Assignment") {
        this->assigner = new Lloyd_assignment<Point>(this->K, this->Grids, this->L, this->k);
    } else if (Assigner == "Inverse Assignment") {
        this->assigner = new Inverse_assignment<Point>(this->K, this->Grids, this->L, this->k);
    } else {
        cerr << "Unknown Assigner";
    }
    cout << '\t' << "Assigner: " << assigner->get_name() << endl;
    /* Updater */
    if (Updater == "Partitioning Around Medoids (PAM)") {
        this->updater = new PAM<Point>(this->K);
    } else if (Updater == "Mean Vector - DTW centroid Curve") {
        this->updater = new MV_DTW<Point>(this->K);
    } else {
        cerr << "Unknown Updater";
    }
    cout << '\t' << "Update: " << updater->get_name() << endl << "-----------------------" << endl;
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset, DistanceDatabase<Point>* db) {
    int convergence = 0;
    int count = 0;
    //for testing = 10
    int max_comps = 10;
    /* initialization */
    this->centroids = initializer->init(dataset);
    /* repeat until no change in clusters */
    cout << endl;
    do {
        cout << "Iteration <" << count << ">: " << endl;
        /* assignment */
        cout << '\t' << "Assigner call ..." << endl;
        this->clusters = assigner->assign(dataset, this->centroids);
        /* update */
        cout << '\t' << "Updater call ..." << endl;
        convergence = updater->update(dataset, this->clusters, this->centroids, db);
        count++;
        /* check convergence */
        if (convergence == 1) break;
        if (count == max_comps)   break;
    }
    while(1);
}

template <class Point>
vector<double> Cluster <Point>::silhouette(vector<vector<Point>>* cluster_data, DistanceDatabase<Point>* db) {
    cout << "Silhouette for cluster statistics ..." << endl;
    double a, b;
    double s = 0;
    vector<double> slt;
    int closest_centroid_id;
    for (int i = 0; i < this->K; i++) {
        s = 0;
        for (int point : *clusters[i]) {
            a = average_distance(point, clusters[i], db);
            closest_centroid_id = find_closest_centroid(centroids[i], db);
            b = average_distance(point, clusters[closest_centroid_id], db);
            s += (b - a)/(double)max(a,b);
        }
        s /= clusters[i]->size();
        slt.push_back(s);
    }
    return slt;
}

template <class Point>
double Cluster <Point>::average_distance(int point, vector<int>* points, DistanceDatabase<Point>* db) {
    double avg = 0.0;
    for (int id : *points) {
        avg += db->get_distance(point, id);
    }
    return avg / (points->size()-1);
}

template <class Point>
int Cluster <Point>::find_closest_centroid(pair<vector<Point>*,int> centroid, DistanceDatabase<Point>* db) {
    double distance = 0.0;
    double min = DBL_MAX;
    int closest_centroid = 0;
    int centroid_id = 0;
    int id;

    for (auto it : this->centroids) {
        id = it.second;
        centroid_id++;
        if ((it.second != -1) && (centroid.second != -1)) {
            if (it.second == centroid.second) continue;
            distance = db->get_distance(centroid.second, id);
        } else {
            distance = dist(centroid.first, it.first, (centroid.first)->size());
        }
        if (distance < min) {
            min = distance;
            closest_centroid = centroid_id-1;
        }
    }
    return closest_centroid;

}

template <class Point>
Cluster <Point>::~Cluster(){
    delete (this->initializer);
}

template class Cluster<double>;
template class Cluster<double*>;