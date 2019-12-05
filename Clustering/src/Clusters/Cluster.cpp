#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
Cluster <Point>::Cluster(int* cluster_conf, string Initializer, string Assigner, string Updater) {

    /* Number of clusters*/
    this->K = cluster_conf[0];
    /* Number of grids*/
    this->Grids = cluster_conf[1];
    /* Number of vector hash tables*/
    this->L = cluster_conf[2];
    /* Number of vector hash functions*/
    this->k = cluster_conf[3];

    cout << "--------Configuration--------" << endl;

    /* Initializer */
    if (Initializer == "Random Selection") {
        this->initializer = new Random_Selection<Point>(this->K);
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(this->K);
    } else {
        cerr << "Unknown Initializer";
    }
    cout << "Initializer: " << initializer->get_name() << endl;

    /* Assigner */
    if (Assigner == "Lloyd's Assignment") {
        this->assigner = new Lloyd_assignment<Point>(this->K, this->Grids, this->L, this->k);
    } else if (Assigner == "Inverse Assignment") {
        this->assigner = new Inverse_assignment<Point>(this->K, this->Grids, this->L, this->k);
    } else {
        cerr << "Unknown Assigner";
    }
    cout  << "Assigner: " << assigner->get_name() << endl;

    /* Updater */
    if (Updater == "Partitioning Around Medoids (PAM)") {
        this->updater = new PAM<Point>(this->K);
    } else if (Updater == "Mean Vector - DTW centroid Curve") {
        this->updater = new MV_DTW<Point>(this->K);
    } else {
        cerr << "Unknown Updater";
    }
    cout << "Updater: " << updater->get_name() << endl;

    cout << "-----------------------------" << endl;
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset, DistanceDatabase<Point>* db) {

    int convergence = 0;
    int count = 0;

    /* Initialization */
    this->centroids = initializer->init(dataset);

    /* Repeat until MAX ITERATIONS or Convergence reached */
    do {
        cout << "\tIteration <" << count+1 << "> " << endl;

        /* Assignment */
        this->clusters = assigner->assign(dataset, &this->centroids, db);

        /* Update */
        convergence = updater->update(dataset, this->clusters, &this->centroids, db);

        count++;

        /* Break Check for convergence */
        if (convergence == 1) {
            /* last assign after the convergence centroids */
            this->clusters = assigner->assign(dataset, &this->centroids, db);
            break;
        }
        /* Break Check for max iterations */
        if (count == MAX_ITERATIONS) {
            /* last assign after the the last update of centroids */
            this->clusters = assigner->assign(dataset, &this->centroids, db);
            break;
        }
    }while(1);

    if( convergence == 1 ){
        cout << " Reached Convergence " << endl;
    }else{
        cout << " Reached MAX ITERATIONS " << endl;
    }
}

template <class Point>
vector<double> Cluster <Point>::silhouette(vector<vector<Point>>* cluster_data, DistanceDatabase<Point>* db) {

    double a, b;
    double s = 0;
    vector<double> slt;
    int closest_centroid_id, cluster_size;
    for (int i = 0; i < this->K; i++) {
        cluster_size = clusters[i]->size();
        s = 0;
        for (int point : *clusters[i]) {
            a = average_distance(point, clusters[i], db);
            closest_centroid_id = find_closest_centroid(centroids[i], db);
            b = average_distance(point, clusters[closest_centroid_id], db);
            s += (b - a)/(double)max(a,b);
        }
        if (cluster_size > 1) {
            s /= cluster_size;
            slt.push_back(s);
        } else {
            slt.push_back(0);
        }
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

        /* for duplicate points in dataset - centroids*/
        if (distance == 0) continue;

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
    delete (this->assigner);

    clear_centroids(&centroids);

    delete (this->updater);

    for (int i = 0; i < this->K; i++)
        delete (this->clusters[i]);
    delete[] this->clusters;
}

template class Cluster<double>;
template class Cluster<double*>;

void clear_centroids(vector<pair<vector<double>*, int>>* centroids) {

    for (int i = 0; i < centroids->size(); i++) {
        if ((*centroids)[i].second == -1)
            delete ((*centroids)[i].first);
    }

    return;
}

void clear_centroids(vector<pair<vector<double*>*, int>>* centroids) {

    for (int i = 0; i < centroids->size(); i++) {
        if ((*centroids)[i].second == -1)
            delete ((*centroids)[i].first);
    }

    return;
}