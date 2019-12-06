#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"
#include "LSH_Functions.h"

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

    /* Database of vector distances in order to be easily reachable when needed and not calculated every time */
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double>* db = new DistanceDatabase<double>();
    db->calculate_distances(&cluster_data);
    cout << "Distances Dictionary completed!" << endl;
    cout << "Starting Program ..." << endl << endl;

    cout << "Average Distance between Neighbors. (Default value for given datasets = 4800)" << endl;
    double w = 4800;
    bool computed_window = false;
    /* pre processed w for the inverse assignment lsh */
    while (!computed_window) {
        char chw;
        cout << "- Press 'W' or 'w' to compute window automatically." << endl << "- Press 'I' or 'i' to insert your value manually below." << endl << "- Press 'D' or 'd' for default value." << endl;
        cin >> chw;
        if (chw == 'w' || chw == 'W') {
            cout << "Computing w ..." << endl;
            w = 4*compute_window(&cluster_data);
            cout << "Proceeding with w = " << w << endl;
            computed_window = true;
        } else if (chw == 'I' || chw == 'i') {
            cout << "Insert w: ";
            cin >> w;
            cout << "Proceeding with w = " << w << endl;
            computed_window = true;
        } else if (chw == 'D' || chw == 'd') {
            cout << "Default value for w" << endl;
            cout << "Proceeding with w = " << w << endl;
            computed_window = true;
        } else {
            cout << "<Unknown command>" << endl;
        }
    }

    ofstream results;
    results.open(results_file);
    for (int i = 0; i < Initializers.size(); i++){
        for(int j = 0; j < Assigners.size(); j++){
            for (int k = 0; k < Updaters.size(); k++){

                string algorithm = "I" + to_string(i+1) + "A" + to_string(j+1) + "U" + to_string(k+1);
                cout << "Algorithm '" << algorithm << "' : IN PROGRESS" << endl;

                /* Clustering Algorithm Execution */
                auto start = chrono::high_resolution_clock::now();
                Cluster <double>* cluster = new Cluster<double>(cluster_config, Initializers[i], Assigners[j], Updaters[k], w);
                cluster->fit(&cluster_data, db);
                auto finish = chrono::high_resolution_clock::now();
                auto elapsed = finish - start;
                double clustering_time = chrono::duration<double>(elapsed).count();

                /* Silhouette Calculation */
                vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);

                /* Centroids and Clusters from Cluster class */
                vector<pair<vector<double>*, int>>* centroids = cluster->get_centroids();
                vector<int>** clusters = cluster->get_clusters();

                /* Clustering Results Output */
                results << "Algorithm: " << algorithm << endl;
                for( int a = 0; a < centroids->size(); a++){
                    results << "CLUSTER-" << a+1 << " {size: " << setw(4) << setfill('0') << clusters[a]->size() << ", centroid: ";
                    if((*centroids)[a].second == -1){
                        results << " { ";
                        for ( int b = 1; b < (*centroids)[a].first->size(); b++){
                            results << setprecision(5) << showpoint << fixed << (*(*centroids)[a].first)[b];
                            if( b+1 != (*centroids)[a].first->size()){
                                results << ", ";
                            }
                        }
                        results << " } " << endl;
                    }else{
                        results << setw(4) << setfill('0') << item_ids[(*centroids)[a].second] << "}" << endl;
                    }
                }
                results << "Clustering Time: " << clustering_time << endl;

                /* Silhouettes Results Output */
                results << "Silhouette: [";
                double silhouette_total = 0;
                for (int i = 0; i < Silhouettes.size(); i++) {
                    results << Silhouettes[i] << ", ";
                    silhouette_total += Silhouettes[i];
                }
                silhouette_total = silhouette_total / Silhouettes.size();
                results << silhouette_total << "]" << endl << endl;

                /* Complete Option Output */
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
                cout << "Algorithm '" << algorithm << "' : COMPLETED SUCCESSFULLY" << endl << endl;
            }
        }
    }
    results.close();
    cout << "Program completed successfully!" << endl;
    delete[] cluster_config;
    delete (db);
    return 0;
}

int Cluster_Curves(string input_file, string config_file, string results_file, int complete){

    int* cluster_config = new int[4];
    vector<vector<double*>> cluster_data;
    vector<string> item_ids;

    /* Read input.dat and cluster.conf and load them in vectors*/
    int error_code = Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids);
    if (error_code == -1) return -1;

    /* Database of curve distances in order to be easily reachable when needed and not calculated every time */
    cout << "Building Distances Dictionary. It might take a while ..." << endl;
    DistanceDatabase<double*>* db = new DistanceDatabase<double*>();
    db->calculate_distances(&cluster_data);
    cout << "Distances Dictionary completed!" << endl;
    cout << "Starting Program ..." << endl << endl;

    cout << "Average Distance between Neighbors. (Default value for given datasets = 300)" << endl;
    double w = 300;
    bool computed_window = false;
    /* pre processed w for the inverse assignment lsh */
    while (!computed_window) {
        char chw;
        cout << "- Press 'I' or 'i' to insert your value manually below." << endl << "- Press 'D' or 'd' for default value." << endl;
        cin >> chw;
        if (chw == 'I' || chw == 'i') {
            cout << "Insert w: ";
            cin >> w;
            cout << "Proceeding with w = " << w << endl;
            computed_window = true;
        } else if (chw == 'D' || chw == 'd') {
            cout << "Default value for w" << endl;
            cout << "Proceeding with w = " << w << endl;
            computed_window = true;
        } else {
            cout << "<Unknown command>" << endl;
        }
    }

    ofstream results;
    results.open(results_file);
    for (int i = 0; i < Initializers.size(); i++){
        for(int j = 0; j < Assigners.size(); j++){
            for (int k = 0; k < Updaters.size(); k++){

                string algorithm = "I" + to_string(i+1) + "A" + to_string(j+1) + "U" + to_string(k+1);
                cout << "Algorithm '" << algorithm << "' : IN PROGRESS" << endl;

                /* Clustering Algorithm Execution */
                auto start = chrono::high_resolution_clock::now();
                Cluster <double*>* cluster = new Cluster<double*>(cluster_config, Initializers[i], Assigners[j], Updaters[k], w);
                cluster->fit(&cluster_data, db);
                auto finish = chrono::high_resolution_clock::now();
                auto elapsed = finish - start;
                double clustering_time = chrono::duration<double>(elapsed).count();

                /* Silhouette Calculation */
                vector<double> Silhouettes = cluster->silhouette(&cluster_data, db);

                /* Centroids and Clusters from Cluster class */
                vector<pair<vector<double*>*, int>>* centroids = cluster->get_centroids();
                vector<int>** clusters = cluster->get_clusters();

                /* Clustering Results Output */
                results << "Algorithm: " << algorithm << endl;
                for( int a = 0; a < centroids->size(); a++){
                    results << "CLUSTER-" << a+1 << " {size: " << setw(4) << setfill('0') << clusters[a]->size() << ", centroid: ";
                    if((*centroids)[a].second == -1){
                        results << " { ";
                        for ( int b = 1; b < (*centroids)[a].first->size(); b++){
                            results << "(" << setprecision(5) << showpoint << fixed << (*(*centroids)[a].first)[b][0] << " , " << setprecision(5) << showpoint << fixed << (*(*centroids)[a].first)[b][1] << ")";
                            if( b+1 != (*centroids)[a].first->size()){
                                results << ", ";
                            }
                        }
                        results << " } " << endl;
                    }else{
                        results << setw(4) << setfill('0') << item_ids[(*centroids)[a].second] << "}" << endl;
                    }
                }
                results << "Clustering Time: " << clustering_time << endl;

                /* Silhouettes Results Output */
                results << "Silhouette: [";
                double silhouette_total = 0;
                for (int i = 0; i < Silhouettes.size(); i++) {
                    results << Silhouettes[i] << ", ";
                    silhouette_total += Silhouettes[i];
                }
                silhouette_total = silhouette_total / Silhouettes.size();
                results << silhouette_total << "]" << endl << endl;

                /* Complete Option Output */
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
                cout << "Algorithm '" << algorithm << "' : COMPLETED SUCCESSFULLY" << endl << endl;
            }
        }
    }
    results.close();
    cout << "Program completed successfully!" << endl;

    /* clear dataset */
    for (auto it : cluster_data){
        for (auto point : it) {
            delete[] point;
        }
    }
    delete[] cluster_config;
    delete (db);
    return 0;
}
