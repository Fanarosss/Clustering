#include "Library.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "BHC.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

int main(int argc, char* argv[]) {

    /* Arguments Check */
    int k = 4, dim = 3, M = 10, probes = 2;             //Default values: k is the number of hi concatenated to form g - dim is number of hypercube's vertices
    double R;                                           //BONUS: Read from Query Input File
    string data_file, search_file, results_file;
    int df = 0, sf = 0, rf = 0;
    char rerun;
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if(arg == "-d"){
            if(i+1 < argc){
                data_file = argv[++i];
                if(access( data_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <input file> -- \n";
                    show_bhc_usage(argv[0]);
                    return -1;
                }
                df = 1;
            }else{
                cout << "-- NO <input file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if (arg == "-q"){
            if(i+1 <= argc){
                search_file = argv[++i];
                if(access( search_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <query file> -- \n";
                    show_bhc_usage(argv[0]);
                    return -1;
                }
                sf = 1;
            }else{
                cout << "-- NO <query file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-o"){
            if(i+1 <= argc){
                results_file = argv[++i];
                rf = 1;
            }else{
                cout << "-- NO <output file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-k"){
            dim = atoi(argv[++i]);
        }else if(arg == "-M") {
            M = atoi(argv[++i]);
        }else if(arg == "-probes"){
            probes = atoi(argv[++i]);
        }else{
            show_bhc_usage(argv[0]);
            return -1;
        }
    }
    /* Loop to rerun program with different files */
    do {

        /* File Check if not given as arguments */
        if (df == 0) {
            cout << "Path to data file (<input file>):" << endl;
            cin >> data_file;
            while (access(data_file.c_str(), F_OK) == -1) {
                cout << "-- Wrong <input file> -- \n";
                cin >> data_file;
            }
        } else df = 0;
        if (sf == 0) {
            cout << "Path to search file (<query file>):" << endl;
            cin >> search_file;
            while (access(search_file.c_str(), F_OK) == -1) {
                cout << "-- Wrong <query file> -- \n";
                cin >> search_file;
            }
        } else sf = 0;
        if (rf == 0) {
            cout << "Path to file of results (<output file>):" << endl;
            cin >> results_file;
        } else rf = 0;

        /* Data & Query vectors */
        vector <vector<int>> dataset;
        vector <vector<int>> searchset;

        /* Read data set and query set and load them in vectors - Also read Radius */
        int error_code = Read_point_files(&dataset, &searchset, &R, data_file, search_file);
        if (error_code == -1) return -1;

        /* Default w*/
        int w = 4 * 1140;
        bool computed_window = false;

        /* Brute Force for actual NNs */
        char bfsearch;
        vector<int> TrueDistances;
        vector<double> TrueTimes;
        cout << endl << "Do you want to run Brute Force? (y/n)\n" ;
        cin >> bfsearch;
        if (bfsearch == 'y' || bfsearch == 'Y') {
            cout << "Starting Exhaustive Search..." << endl;
            brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes);
            cout << "Exhaustive Search completed." << endl;
        }else{
            string brute_force_file;
            cout << "Path to brute force file:" << endl;
            cin >> brute_force_file;
            while (access(brute_force_file.c_str(), F_OK) == -1) {
                cout << "-- Wrong brute force file -- \n";
                cin >> brute_force_file;
            }
            read_vectors_brute_force_file(brute_force_file, &TrueDistances, &TrueTimes);
        }

        /* Arrays for results */
        int *min_distance = new int[searchset.size()];
        int *nearest_neighbor = new int[searchset.size()];
        double *time = new double[searchset.size()];

        /* Initialize arrays */
        for (int i = 0; i < searchset.size(); i++) {
            min_distance[i] = INT_MAX;
            nearest_neighbor[i] = -1;
            time[i] = 0;
        }

        /* Compute w */
        while (!computed_window) {
            char chw;
            cout << "- Press 'W' or 'w' to compute window automatically." << endl << "- Press 'I' or 'i' to insert your value manually below." << endl << "- Press 'D' or 'd' for default value." << endl;
            cin >> chw;
            if (chw == 'w' || chw == 'W') {
                cout << "Computing w ..." << endl;
                w = 4*compute_window(&dataset);
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

        /* Vector for R-Neighbors (BONUS) */
        vector<vector<int>> R_neighbors;

        /* ---- CALL BHC ---- */
        BHC(&dataset, &searchset, k, dim, M, probes, w, R, &R_neighbors, &min_distance, &time, &nearest_neighbor);

        /* Variables for effectiveness and time*/
        double max_af = 0.0;
        double average_af = 0.0;
        double curr_fraction = 0.0;
        double average_time = 0.0;

        /* Results for every query */
        for (int q = 0; q < searchset.size(); q++) {
            curr_fraction = (double) min_distance[q] / TrueDistances[q];
            if (curr_fraction > max_af) max_af = curr_fraction;
            average_af += curr_fraction;
            average_time += time[q];
        }
        average_af = average_af / searchset.size();
        average_time = average_time / searchset.size();

        /* Print used variables and results */
        cout << "Variables used: | k = " << k << " | dim = " << dim << " | M = " << M << " | probes = " << probes
             << endl;
        cout << "MAX Approximation Fraction (BHC Distance / True Distance) = " << max_af << endl;
        cout << "Average Approximation Fraction (BHC Distance / True Distance) = " << average_af << endl;
        cout << "Average Time of BHC Distance Computation = " << setprecision(9) << showpoint << fixed << average_time
             << endl;

        /* Write results in file */
        ofstream neighbors_file;
        neighbors_file.open("./output/" + results_file);
        for (int i = 0; i < searchset.size(); i++) {
            neighbors_file << "Query: " << i + 1 << endl;
            neighbors_file << "Nearest Neighbor: " << nearest_neighbor[i] + 1 << endl;
            neighbors_file << "distanceLSH: " << min_distance[i] << endl;
            neighbors_file << "distanceTrue: " << TrueDistances[i] << endl;
            neighbors_file << "tLSH: " << setprecision(9) << showpoint << fixed << time[i] << endl;
            neighbors_file << "tTrue: " << setprecision(9) << showpoint << fixed << TrueTimes[i] << endl;
            if (R != 0) {
                neighbors_file << "R-near neighbors: " << endl;
                if (R_neighbors[i].size() != 0) {
                    for (int j = 0; j < R_neighbors[i].size(); j++) {
                        neighbors_file << R_neighbors[i][j] + 1 << endl;
                    }
                } else {
                    neighbors_file << "No R-near neighbors available" << endl;
                }
            }
            neighbors_file << endl;
        }
        neighbors_file.close();

        /* Clean Memory */
        delete[] min_distance;
        delete[] nearest_neighbor;
        delete[] time;

        /* Clear underlying memory of vectors for next iteration */
        vector<vector<int>>().swap(R_neighbors);
        vector <vector<int>>().swap(dataset);
        vector <vector<int>>().swap(searchset);
        vector<int>().swap(TrueDistances);
        vector<double>().swap(TrueTimes);

        /* Option to rerun program with different files */
        cout<<"\nDo you want to run this program again? (y/n)\n";
        cin>>rerun;
    }while (rerun == 'y' || rerun == 'Y');

    return 0;
}
