#include "Library.h"
#include "Cluster.h"

using namespace std;

template <class Point>
Cluster <Point>::Cluster(int K, string Initializer, string Assigner, string Updater) {
    /* k for k means */
    this->K = K;
    /* Initializer */
    if (Initializer == "Random Selection") {
        this->initializer = new Random_Selection<Point>(K);
        cout << initializer->get_name() << endl;
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(K);
        cout << initializer->get_name() << endl;
    } else {
        cerr << "Unknown Initializer";
    }
    /* Assigner */
    if (Assigner == "Lloyd's Assignment") {
        this->assigner = new Lloyd_assignment<Point>(K);
        cout << assigner->get_name() << endl;
    } else if (Assigner == "Inverse Assignment") {
        this->assigner = new Inverse_assignment<Point>(K);
        cout << assigner->get_name() << endl;
    } else {
        cerr << "Unknown Assigner";
    }
    /* Updater */
    if (Updater == "Partitioning Around Medoids (PAM)") {
        this->updater = new PAM<Point>(K);
        cout << updater->get_name() << endl;
    } else if (Updater == "Mean Vector - DTW centroid Curve") {
        this->updater = new MV_DTW<Point>(K);
        cout << updater->get_name() << endl;
    } else {
        cerr << "Unknown Updater";
    }
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset) {
    /* initialization */
    cout << '\t' << "Initializer call ..." << endl;
    this->centroids = initializer->init(dataset);
    /* assignment */
    cout << '\t' << "Assigner call ..." << endl;
    this->clusters = assigner->assign(dataset, this->centroids);
    /* update */
    cout << '\t' << "Updater call ..." << endl;
    updater->update(dataset, this->clusters, this->centroids);
}

template <class Point>
Cluster <Point>::~Cluster(){
    delete (this->initializer);
}

template class Cluster<int>;
template class Cluster<double*>;