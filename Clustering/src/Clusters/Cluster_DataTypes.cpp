#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"

using namespace std;

int Cluster_Vectors(string input_file, string config_file){
    int* cluster_config = new int[4];
    vector<vector<double>> cluster_data;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file);
    if (error_code == -1) return -1;
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double>* db = new DistanceDatabase<double>();
    db->calculate_distances(&cluster_data);
    cout << "Clustering vectors..." << endl;
    string initializer = "K-Means++";
    string assigner = "Inverse Assignment";
    string updater = "Mean Vector - DTW centroid Curve";
    Cluster <double>* cluster = new Cluster<double>(cluster_config, initializer, assigner, updater);
    cluster->fit(&cluster_data, db);
    vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);
    /* printing Silhouettes */
    for (int i = 0; i < Silhouettes.size(); i++) {
        cout << '\t' << "Cluster <" << i << "> : " << Silhouettes[i] << endl;
    }
    delete (cluster);
    delete[] cluster_config;
    return 0;
}

int Cluster_Curves(string input_file, string config_file){
    int* cluster_config = new int[4];
    vector<vector<double*>> cluster_data;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file);
    if (error_code == -1) return -1;
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double*>* db = new DistanceDatabase<double*>();
    db->calculate_distances(&cluster_data);
    cout << "Clustering vectors..." << endl;
    string initializer = "K-Means++";
    string assigner = "Lloyd's Assignment";
    string updater = "Mean Vector - DTW centroid Curve";
    Cluster <double*>* cluster = new Cluster<double*>(cluster_config, initializer, assigner, updater);
    cluster->fit(&cluster_data, db);
    vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);
    /* printing Silhouettes */
    for (int i = 0; i < Silhouettes.size(); i++) {
        cout << '\t' << "Cluster <" << i << "> : " << Silhouettes[i] << endl;
    }
    delete (cluster);
    delete[] cluster_config;
    return 0;
}
