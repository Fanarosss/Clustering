#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"

using namespace std;

map<int,string> Initializers = {{0 , "Random Selection"} , {1 , "K-Means++"}};
map<int,string> Assigners = {{0 , "Lloyd's Assignment"} , {1 , "Inverse Assignment"}};
map<int,string> Updaters = {{0 , "Partitioning Around Medoids (PAM)"} , {1 , "Mean Vector - DTW centroid Curve"}};

int Cluster_Vectors(string input_file, string config_file, string results_file, int complete){
    int* cluster_config = new int[4];
    vector<vector<double>> cluster_data;
    vector<string> item_ids;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids);
    if (error_code == -1) return -1;
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double>* db = new DistanceDatabase<double>();
    db->calculate_distances(&cluster_data);
    cout << "Distances Dictionary completed!" << endl;
    cout << " Starting Program ..." << endl << endl;

    ofstream results;
    results.open(results_file);
    for (int i = 0; i < Initializers.size(); i++){
        for(int j = 0; j < Assigners.size(); j++){
            for (int k = 0; k < Updaters.size(); k++){
                string algorithm = "I" + to_string(i+1) + "A" + to_string(j+1) + "U" + to_string(k+1);
                cout << "Algorithm '" << algorithm << "' : IN PROGRESS" << endl;
                auto start = chrono::high_resolution_clock::now();
                Cluster <double>* cluster = new Cluster<double>(cluster_config, Initializers[i], Assigners[j], Updaters[k]);
                cluster->fit(&cluster_data, db);
                auto finish = chrono::high_resolution_clock::now();
                auto elapsed = finish - start;
                double clustering_time = chrono::duration<double>(elapsed).count();
                vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);

                vector<pair<vector<double>*, int>>* centroids = cluster->get_centroids();
                vector<int>** clusters = cluster->get_clusters();

                results << "Algorithm: " << algorithm << endl;
                for( int a = 0; a < centroids->size(); a++){
                    results << "CLUSTER-" << a+1 << " {size: " << clusters[a]->size() << ", centroid: ";
                    if((*centroids)[a].second == -1){
                        results << " { ";
                        for ( int b = 1; b < (*centroids)[a].first->size(); b++){
                            results << (*(*centroids)[a].first)[b];
                            if( b+1 != (*centroids)[a].first->size()){
                                results << ", ";
                            }
                        }
                        results << " } " << endl;
                    }else{
                        results << (*centroids)[a].second << "}" << endl;
                    }
                }
                results << "clustering time: " << clustering_time << endl;
                /* printing Silhouettes */
                results << "Silhouette: [";
                double silhouette_total = 0;
                for (int i = 0; i < Silhouettes.size(); i++) {
                    results << Silhouettes[i] << ", ";
                    silhouette_total += Silhouettes[i];
                }
                silhouette_total = silhouette_total / Silhouettes.size();
                results << silhouette_total << "]" << endl << endl;
                if( complete == 1 ){
                    for( int a = 0; a < centroids->size(); a++) {
                        results << "CLUSTER-" << a+1 << " { ";
                        for ( int b = 0; b < clusters[a]->size(); b++){
                            results << item_ids[(*clusters[a])[b]];
                            if( b+1 != clusters[a]->size()){
                                results << ", ";
                            }
                        }
                        results << " } " << endl;
                    }
                    results << endl;
                }


                delete (cluster);
                cout << "Algorithm '" << algorithm << "' : COMPLETED SUCCESSFULLY" << endl;
            }
        }
    }
    results.close();
    cout << endl << "Program completed successfully!" << endl;
    delete[] cluster_config;
    delete (db);
    return 0;

//    cout << "Clustering vectors..." << endl;
//    string initializer = "K-Means++";
//    string assigner = "Lloyd's Assignment";
//    string updater = "Mean Vector - DTW centroid Curve";
//    Cluster <double>* cluster = new Cluster<double>(cluster_config, initializer, assigner, updater);
//    cluster->fit(&cluster_data, db);
//    vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);
//    /* printing Silhouettes */
//    for (int i = 0; i < Silhouettes.size(); i++) {
//        cout << '\t' << "Cluster <" << i << "> : " << Silhouettes[i] << endl;
//    }
//    delete (cluster);
//    delete[] cluster_config;
//    delete (db);
//    return 0;
}

int Cluster_Curves(string input_file, string config_file, string results_file, int complete){
    int* cluster_config = new int[4];
    vector<vector<double*>> cluster_data;
    vector<string> item_ids;
    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids);
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
