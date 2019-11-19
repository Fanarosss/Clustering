#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"

using namespace std;

int Cluster_Vectors(string input_file, string config_file){
    int* cluster_config = new int[4];
    vector<vector<int>> cluster_data;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file);
    if (error_code == -1) return -1;

    cout << "Clustering vectors..." << endl;
    string initializer = "K-Means++";
    string assigner = "Lloyd's Assignment";
    string updater = "Partitioning Around Medoids (PAM)";
    Cluster <int>* cluster = new Cluster<int>(5, initializer, assigner, updater);
    cluster->fit(&cluster_data);
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

    cout << "Clustering vectors..." << endl;
    string initializer = "K-Means++";
    string assigner = "Lloyd's Assignment";
    string updater = "Partitioning Around Medoids (PAM)";
    Cluster <double*>* cluster = new Cluster<double*>(5, initializer, assigner, updater);
    // cluster->fit(&cluster_data);
    delete (cluster);
    delete[] cluster_config;
    return 0;
}
