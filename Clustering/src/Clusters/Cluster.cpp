#include "Library.h"
#include "Cluster.h"

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
        this->initializer = new Random_Selection<Point>(K);
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(K);
    } else {
        cerr << "Unknown Initializer";
    }
    cout << '\t' << "Initializer: " << initializer->get_name() << endl;
    /* Assigner */
    if (Assigner == "Lloyd's Assignment") {
        this->assigner = new Lloyd_assignment<Point>(K);
    } else if (Assigner == "Inverse Assignment") {
        this->assigner = new Inverse_assignment<Point>(K);
    } else {
        cerr << "Unknown Assigner";
    }
    cout << '\t' << "Assigner: " << assigner->get_name() << endl;
    /* Updater */
    if (Updater == "Partitioning Around Medoids (PAM)") {
        this->updater = new PAM<Point>(K);
    } else if (Updater == "Mean Vector - DTW centroid Curve") {
        this->updater = new MV_DTW<Point>(K);
    } else {
        cerr << "Unknown Updater";
    }
    cout << '\t' << "Update: " << updater->get_name() << endl << "-----------------------" << endl;
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset) {
    int convergence = 0;
    int count = 0;
    //for testing = 1
    int max_comps = 1;
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
        convergence = updater->update(dataset, this->clusters, this->centroids);
        count++;
        /* check convergence */
        if (convergence == 1) break;
        if (count == max_comps)   break;
    }
    while(1);
}

template <class Point>
Cluster <Point>::~Cluster(){
    delete (this->initializer);
}

template class Cluster<int>;
template class Cluster<double*>;