#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"

using namespace std;

map<int,string> Initializers = {{0 , "Random Selection"} , {1 , "K-Means++"}};
map<int,string> Assigners = {{0 , "Lloyd's Assignment"} , {1 , "Inverse Assignment"}};
map<int,string> Updaters = {{0 , "Partitioning Around Medoids (PAM)"} , {1 , "Mean Vector - DTW centroid Curve"}};

int Cluster_Vectors(string input_file, string config_file, string results_file){
    int* cluster_config = new int[4];
    vector<vector<double>> cluster_data;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file);
    if (error_code == -1) return -1;
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double>* db = new DistanceDatabase<double>();
    db->calculate_distances(&cluster_data);
    cout << "Distances Dictionary completed! Starting Program ..." << endl << endl;

//    for (int i = 0; i < Initializers.size(); i++){
//        for(int j = 0; j < Assigners.size(); j++){
//            for (int k = 0; k < Updaters.size(); k++){
//                string algorithm = "I" + to_string(i+1) + "A" + to_string(j+1) + "U" + to_string(k+1);
//                cout << "Algorithm '" << algorithm << "' : IN PROGRESS" << endl;
//                Cluster <double>* cluster = new Cluster<double>(cluster_config, Initializers[i], Assigners[j], Updaters[k]);
//                cluster->fit(&cluster_data, db);
//                vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);
//                /* printing Silhouettes */
//                for (int i = 0; i < Silhouettes.size(); i++) {
//                    cout << '\t' << "Cluster <" << i << "> : " << Silhouettes[i] << endl;
//                }
//                delete (cluster);
//                cout << "Algorithm '" << algorithm << "' : COMPLETED SUCCESSFULLY" << endl;
//            }
//        }
//    }
//    delete[] cluster_config;
//    delete (db);
//    return 0;

    cout << "Clustering vectors..." << endl;
    string initializer = "K-Means++";
    string assigner = "Lloyd's Assignment";
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
    delete (db);
    return 0;
}

int Cluster_Curves(string input_file, string config_file, string results_file){
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
    delete (db);
    for (auto curve : cluster_data){
        for (auto point : curve){
            delete[] point;
        }
    }
    return 0;
}
