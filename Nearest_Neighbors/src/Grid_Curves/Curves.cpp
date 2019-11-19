#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"
#include "BHC.h"
#include "LSH.h"
#include "LSH_Functions.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* variable declaration | L_grid = 4 default value */
    int L_grid = 4;
    double R = 0;
    bool computed_window = false;
    /* default 2D curves */
    int d = 2;

    /* default values for LSH and Hypercube */
#ifdef _LSH_
    int k_vec = 1, L_vec = 1;
#endif

#ifdef _BHC_
    int k_hypercube = 1, M = 10, probes = 10, dim = 3;
#endif

    /* read arguments */
    string data_file, search_file, results_file;
    int df = 0, sf = 0, rf = 0, input_flag = 0;
    char rerun;
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if(arg == "-d"){
            if(i+1 < argc){
                data_file = argv[++i];
                if(access( data_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <input file> -- \n";
                    input_flag = 1;
                    break;
                }
                df = 1;
            }else{
                cout << "-- NO <input file> -- \n";
                input_flag = 1;
                break;
            }
        }else if (arg == "-q"){
            if(i+1 <= argc){
                search_file = argv[++i];
                if(access( search_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <query file> -- \n";
                    input_flag = 1;
                    break;
                }
                sf = 1;
            }else{
                cout << "-- NO <query file> -- \n";
                input_flag = 1;
                break;
            }
        }else if(arg == "-o"){
            if(i+1 <= argc){
                results_file = argv[++i];
                rf = 1;
            }else{
                cout << "-- NO <output file> -- \n";
                input_flag = 1;
                break;
            }
        }else if(arg == "-L_grid"){
            L_grid = atoi(argv[++i]);
        }

#ifdef _LSH_
        else if(arg == "-k_vec"){
            k_vec = atoi(argv[++i]);
        }
#endif

#ifdef _BHC_
        /* k_hypercube(=d') as stated in the discussions */
        else if(arg == "-k_hypercube"){
            dim = atoi(argv[++i]);
        }
        else if(arg == "-M"){
            M = atoi(argv[++i]);
        }
        else if(arg == "-probes"){
            probes = atoi(argv[++i]);
        }
#endif

        else{
            input_flag = 1;
            break;
        }
    }

    /* --Help: missing data -> print usage */
    if (input_flag == 1) {
#ifdef _LSH_
        show_grid_lsh_usage(argv[0]);
#endif
#ifdef _BHC_
        show_grid_bhc_usage(argv[0]);
#endif
        return -1;
    }
    /* end of reading data */

    /* verify data */
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
    /* end of verifying data */

    /* start of curve hashing */
    int m1, m2, error_code, min;
    double delta;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    vector<vector<double*>> searchset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files(&dataset, &searchset, data_file, search_file);
    if (error_code == -1){
        cout << "Exiting ..." << endl;
        return -1;
    }

    /* dataset sizes */
    m1 = dataset.size();
    m2 = searchset.size();
    /* calculate delta */
    min = (m1 < m2) ? m1 : m2;
    /* delta = very small value to slice the grid in relevant with our dataset chunks */
    delta = 0.00006;

    double w = 120; /* default value that seems to work quite well */
    double max_element = 0.0;
    int max_points = 0, elements = 0;

    /* do brute force to find actual NNs */
    char bfsearch;
    vector<double> TrueDistances;
    vector<double> TrueTimes;
    vector<int> TrueNeighbors;
    cout << endl << "Do you want to run Brute Force? (y/n)\n" ;
    cin >> bfsearch;
    if (bfsearch == 'y' || bfsearch == 'Y') {
        cout << "Exhaustive Search using DTW. It might take a while ..." << endl;
        curves_brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes, &TrueNeighbors);
        cout << "Found exact neighbors." << endl;
    }else{
        string brute_force_file;
        cout << "Path to brute force file:" << endl;
        cin >> brute_force_file;
        while (access(brute_force_file.c_str(), F_OK) == -1) {
            cout << "-- Wrong brute force file -- \n";
            cin >> brute_force_file;
        }
        read_curves_brute_force_file(brute_force_file, &TrueDistances, &TrueTimes, &TrueNeighbors);
    }

    /* orthogonal grid of size d */
    vector<double> orthogonal_grid;
    /* vector for hashed curves */
    vector<vector<vector<double>>> hashed_curves;
    vector<vector<double>> temp_hash;
    /* vectored curves */
    vector<vector<double>> data_vectored_curves;
    vector<vector<double>> search_vectored_curves;
    /* temp */
    vector<double> curve;
    /* RESULTS storing */
    double *min_distance;
    int *nearest_neighbor;
    double *time;
    double average_time = 0.0;
    /* first hashing */
    vector<int*> hashed_neighbors;
    /* bonus r radius */
    vector<vector<int>> R_neighbors;

    cout << "Running with default delta = " << delta << endl << "Proceeding to preprocessing and hashing ..." << endl;
    /*  ----------- Loop this L times and then dtw on those L nn sets -------- */
    for (int i = 0; i < L_grid; i++) {
        /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
        shift_grid(&orthogonal_grid, delta, d);

        /* ------------------ DATA SET hashing ----------------- */
        /* hash all dataset curves */
        for (int i = 0; i < dataset.size(); i++) {
            hash_curve(&temp_hash, &dataset[i], &orthogonal_grid, delta, d);
            hashed_curves.push_back(temp_hash);
            /* clean temp hash */
            vector<vector<double>>().swap(temp_hash);
        }

        /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
         * the equivalent points in our new grid, that the polygonal curve projects */
        elements = 0;
        max_points = 0;
        max_element = 0.0;
        /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
        for (int i = 0; i < hashed_curves.size(); i++) {
            for (int j = 0; j < hashed_curves[i].size(); j++) {
                curve.push_back(hashed_curves[i][j][0]);
                /* find max element */
                /* j != 0 because at index 0 there is the id and length of curve */
                if (j != 0) {
                    curve.push_back(hashed_curves[i][j][1]);
                    if (hashed_curves[i][j][0] > max_element) {
                        max_element = hashed_curves[i][j][0];
                    }
                    if (hashed_curves[i][j][1] > max_element) {
                        max_element = hashed_curves[i][j][1];
                    }
                }
                elements++;
            }

            if (elements > max_points) {
                max_points = elements;
            }
            elements = 0;
            data_vectored_curves.push_back(curve);
            vector<double>().swap(curve);
        }

        /* ------------------ SEARCH SET hashing ----------------- */
        /* clean hashed curves */
        vector<vector<vector<double>>>().swap(hashed_curves);
        /* end of cleaning */

        /* hash all curves */
        for (int i = 0; i < searchset.size(); i++) {
            hash_curve(&temp_hash, &searchset[i], &orthogonal_grid, delta, d);
            hashed_curves.push_back(temp_hash);
            /* clean temp hash */
            vector<vector<double>>().swap(temp_hash);
        }

        /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
         * the equivalent points in our new grid, that the polygonal curve projects */

        /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
        elements = 0;
        for (int i = 0; i < hashed_curves.size(); i++) {
            for (int j = 0; j < hashed_curves[i].size(); j++) {
                /* we always push back the first element - even the id on (0,0)
                 * we won't push back the dimensions of a curve in (0,1) */
                curve.push_back(hashed_curves[i][j][0]);
                /* j != 0 because at index 0 there is the id and length of curve */
                if (j != 0) {
                    /* we dont push back the length size for the lsh */
                    curve.push_back(hashed_curves[i][j][1]);
                    if (hashed_curves[i][j][0] > max_element) {
                        max_element = hashed_curves[i][j][0];
                    }
                    if (hashed_curves[i][j][1] > max_element) {
                        max_element = hashed_curves[i][j][1];
                    }
                }
                elements++;
            }
            if (elements > max_points) {
                max_points = elements;
            }
            elements = 0;
            search_vectored_curves.push_back(curve);
            vector<double>().swap(curve);
        }

        /* clear no longer used vectors for memory optimzations */
        vector<double>().swap(orthogonal_grid);
        /* clean hashed curves */
        vector<vector<vector<double>>>().swap(hashed_curves);
        /* end of cleaning */

        /* ----------------- PADDING for both sets ------------ */
        /* pad special number > max coord */
        for (int i = 0; i < data_vectored_curves.size(); i++) {
            while (data_vectored_curves[i].size() < max_points) {
                data_vectored_curves[i].push_back(2 * max_element);
            }
        }
        for (int i = 0; i < search_vectored_curves.size(); i++) {
            while (search_vectored_curves[i].size() < max_points) {
                search_vectored_curves[i].push_back(2 * max_element);
            }
        }
        /* ----------- end of padding ------------ */

        /* --- initializations for LSH --- */
        min_distance = new double [searchset.size()];
        nearest_neighbor = new int [searchset.size()];
        time = new double [searchset.size()];

        for (int i = 0; i < searchset.size(); i++) {
            min_distance[i] = -1;
            nearest_neighbor[i] = -1;
            time[i] = 0;
        }

        while (!computed_window) {
            char chw;
            cout << "- Press 'W' or 'w' to compute window automatically." << endl << "- Press 'I' or 'i' to insert your value manually below." << endl << "- Press 'D' or 'd' for default value." << endl;
            cin >> chw;
            if (chw == 'w' || chw == 'W') {
                cout << "Computing w ..." << endl;
                w = 4*compute_window(&data_vectored_curves);
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

#ifdef _BHC_
        /* ---------------- Hashing them again with Hypercube ------------------ */
        BHC(&data_vectored_curves, &search_vectored_curves, k_hypercube, dim, M, probes, w, R, &R_neighbors, &min_distance, &time, &nearest_neighbor);
#endif

#ifdef _LSH_
        /* ---------------- Hashing them again with LSH ------------------ */
         LSH <double>* model = new LSH <double> (k_vec, L_vec, w);
        model->fit(&data_vectored_curves);
        model->evaluate(&search_vectored_curves, R, &R_neighbors, &min_distance, &time, &nearest_neighbor);
        delete (model);
#endif

        /* store results for all iterations of hashing */
        hashed_neighbors.push_back(nearest_neighbor);
        /* compute average_time */
        double sum_times = 0.0;
        for (int s = 0; s < searchset.size(); s++)
            sum_times += time[s];
        average_time += sum_times / searchset.size();
        /* clean vectors for next iteration */
        vector<vector<double>>().swap(data_vectored_curves);
        vector<vector<double>>().swap(search_vectored_curves);
        vector<vector<int>>().swap(R_neighbors);
        /* clean pointers */
        delete[] min_distance;
        delete[] time;
    }

    /* min distance vector is for the lsh hashed data, we will do DTW now on the real curves to find the
     * true distance between the approximate nearest neighbors found by lsh*/
    double distance = 0.0;
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;

    /* initialize arrays to keep the results */
    min_distance = new double[searchset.size()];
    nearest_neighbor = new int[searchset.size()];
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = -1;
        nearest_neighbor[i] = -1;
    }

    /* search on the results fo the LSH\Hypercube with DTW
     * for the closest of the L nearest neighbors */
    for (int i = 0; i < searchset.size(); i++) {
       for (int j = 0; j < L_grid; j++) {
           if (hashed_neighbors[j][i] == -1) continue;
           auto start = chrono::high_resolution_clock::now();
           distance = DTW(&searchset[i], &dataset[hashed_neighbors[j][i]]);
           if (j == 0) {
               min_distance[i] = distance;
               nearest_neighbor[i] = hashed_neighbors[j][i];
           } else if (distance < min_distance[i]) {
               min_distance[i] = distance;
               nearest_neighbor[i] = hashed_neighbors[j][i];
           }
           auto finish = chrono::high_resolution_clock::now();
           auto elapsed = finish - start;
           average_time += chrono::duration<double>(elapsed).count() / searchset.size();
       }
    }

    /* computing the statistics */
    for (int q = 0; q < searchset.size(); q++) {
        if (min_distance[q] == -1 || TrueDistances[q] == 0) continue;
        curr_fraction = (double) min_distance[q] / TrueDistances[q];
        if (curr_fraction > max_af) max_af = curr_fraction;
        average_af += curr_fraction / searchset.size();
    }

    /* --- RESULTS --- */
    average_time = average_time / searchset.size();
    cout << "MAX Approximation Fraction (Grid Curves Distance / True Distance) = " << max_af << endl;
    cout << "Average Approximation Fraction (Grid Curves Distance / True Distance) = " << average_af << endl;
    cout << "Average Time for Nearest Neighbor = " << average_time << endl;

    /* open file to write results */
    string Method;
    #ifdef _LSH_
    Method = "LSH";
    #endif
    #ifdef _BHC_
    Method = "HyperCube";
    #endif
    ofstream neighbors_file;
    /* open file to dump all query results */
    neighbors_file.open("./output/" + results_file);
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Query: " << i + 1 << endl;
        neighbors_file << "Method: LSH" << endl;
        neighbors_file << "HashFunction: " << Method << endl;
        if (nearest_neighbor[i] != -1) {
            neighbors_file << "Found Nearest Neighbor: " << nearest_neighbor[i] + 1 << endl;
        }else{
            neighbors_file << "Found Nearest Neighbor: Fail" << endl;
        }
        neighbors_file << "True Nearest Neighbor: " << TrueNeighbors[i]<< endl;
        if (nearest_neighbor[i] != -1) {
            neighbors_file << "distanceFound: " << min_distance[i] << endl;
        }else{
            neighbors_file << "distanceFound: None" << endl;
        }
        neighbors_file << "distanceTrue: " << TrueDistances[i] << endl;
        neighbors_file << endl;
    }
    neighbors_file.close();

    cout << "Wrote results on " << results_file << endl;

    /* clean remaining used memory */
    delete[] min_distance;
    delete[] nearest_neighbor;
    for (int i = 0; i < hashed_neighbors.size(); i++)
            delete[] hashed_neighbors[i];
    for (int i = 0; i < dataset.size(); i++)
        for (int j = 0; j < dataset[i].size(); j++)
            delete[] dataset[i][j];
    for ( int i = 0; i  < searchset.size(); i++)
        for (int j = 0; j < searchset[i].size(); j++)
            delete[] searchset[i][j];
    vector<int*>().swap(hashed_neighbors);
    /* end of cleaning */

    return 0;
}
