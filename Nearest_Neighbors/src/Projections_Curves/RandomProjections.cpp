#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"
#include "LSH.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* check arguments */
    if (argc != 3) {
        cout << "We need input_file and query file" << endl;
        return -1;
    }

    /* variable declaration | k = 4 default value */
    int k = 1, L = 5, d = 2, L_vec = 1;
    double epsilon = 0.5;
    int m1, m2, error_code, min;
    double delta;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    vector<vector<double*>> searchset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files_max_dim(&dataset, &searchset, argv[1], argv[2], 5.00);
    if (error_code == -1) return -1;

    /* dataset sizes */
    m1 = dataset.size();
    m2 = searchset.size();
    /* calculate delta */
    min = (m1 < m2) ? m1 : m2;
    //delta = 4*d*min - 1;
    delta = 0.12;

    /* ------------------------ HASHING with RANDOM PROJECTIONS ----------------------- */

    /* Let M be the max length of all input and possible query curves */
    int M = 0;
    for (int i = 0; i < dataset.size(); i++) {
        if (dataset[i][0][1] > M) {
            M = dataset[i][0][1];
        }
    }
    for (int i = 0; i < searchset.size(); i++) {
        if (searchset[i][0][1] > M) {
            M = searchset[i][0][1];
        }
    }

    /* do brute force to find actual NNs */
    char bfsearch;
    vector<double> TrueDistances;
    vector<double> TrueTimes;
    vector<int> TrueNeighbors;
    cout << "Do you want to run Brute Force for extra statistics at the end?" << endl << "(Press y or Y + enter to run, else press any other character)." << endl;
    cin >> bfsearch;
    if (bfsearch == 'y' || bfsearch == 'Y') {
        cout << "Exhaustive Search using DTW. It might take a while ..." << endl;
        curves_brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes, &TrueNeighbors);
        cout << "Found exact neighbors." << endl;
    }

    /* ------------------ RELEVANT traversals ----------------- */

    cout << "Creating MxM array where M is: " << M << endl;
    /* Find relevant traversals, all the pairs (Ui,Vi)
     * and we keep the indexes for every pair of curve */
    vector<vector<vector<vector<double>>>>** TraversalsTable = new vector<vector<vector<vector<double>>>>* [M];
    for (int i = 0; i < M; i++) {
        TraversalsTable[i] = new vector<vector<vector<vector<double>>>> [M];
    }

    vector<vector<vector<int>>> traversals;
    vector<vector<vector<vector<double>>>> traversals_to_coordinates;
    vector<vector<vector<double>>> traverse_with_coords;
    vector<vector<double>> pair_coords;
    vector<double> coords;
    int index_x;
    int index_y;
    /* for every pair of curves */
    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < searchset.size(); j++) {
            /* find all relevant traversals */
            Relevant_Traversals(&traversals, dataset[i].size(), searchset[j].size());
            /* Convert the vector of pair with indexes, to a vector with coordinates */
            /* for all traversals of size */
            for (int k = 0; k < traversals.size(); k++) {
                /* for every pair */
                vector<double> pair_ids;
                pair_ids.push_back((double)i);
                pair_ids.push_back((double)j);
                pair_coords.push_back(pair_ids);
                traverse_with_coords.push_back(pair_coords);
                vector<vector<double>>().swap(pair_coords);
                for (int t = 0; t < traversals[k].size(); t++) {
                    /* U */
                    index_x = traversals[k][t][0];
                    index_y = traversals[k][t][1];
                    coords.push_back(dataset[i][index_x][0]);
                    coords.push_back(dataset[i][index_x][1]);
                    pair_coords.push_back(coords);
                    vector<double>().swap(coords);
                    /* V */
                    coords.push_back(searchset[j][index_y][0]);
                    coords.push_back(searchset[j][index_y][1]);
                    pair_coords.push_back(coords);
                    vector<double>().swap(coords);
                    /* push pair */
                    traverse_with_coords.push_back(pair_coords);
                    vector<vector<double>>().swap(pair_coords);
                }
                /* push back the traversal */
                TraversalsTable[(int)dataset[i][0][1]-1][(int)searchset[j][0][1]-1].push_back(traverse_with_coords);
                /* clear vector */
                vector<vector<vector<double>>>().swap(traverse_with_coords);
            }
            vector<vector<vector<int>>>().swap(traversals);
        }
    }

    cout << "Array MxM is ready" << endl;

    /* Let matrix G have real, independent, normally distributed elements ∼N(0,1)
     * and dimensions K × d, where K = −d*log(epsilon) / epsilon2 */
    double K = -d*log2(epsilon)/pow(epsilon,2);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    /* uniformly values in [0,1) */
    vector<vector<double>> G;
    vector<double> G_temp;
    normal_distribution<double> distribution (0, 1);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < d; j++) {
            G_temp.push_back(distribution(generator));
        }
        G.push_back(G_temp);
        G_temp.clear();
        G_temp.shrink_to_fit();
    }
    cout << "Vector G is ready" << endl;

    /* -- LSH structures -- */
    double w = 400;
    double R = 0;
    /* results */
    double *min_distance;
    int *nearest_neighbor;
    double *time;
    /* results for bonus */
    vector<vector<int>> R_neighbors;
    /* ids of traversals */
    vector<vector<vector<int>>> traversal_neighbors;

    /* We project the Ui’s to vector x = [G·U1|···|G·Uu] ∈ R^uK
     * by multiplying G with every point vector, then concatenating. */
    vector<double> termU;
    vector<double> termV;
    vector<vector<double>> Vectored_Traversals_X;
    vector<vector<double>> Vectored_Traversals_Y;
    vector<double> Map_TraversalX_to_Curve;
    vector<double> Map_TraversalY_to_Curve;

    double Unew, Vnew;
    /* traversal counter */
    double traversal_counter = 0;
    /* for 1st dim of MxM table */
    for (int i = 0; i < M; i++) {
        traversal_counter = 0.0;
        /* for 2nd dim of MxM  table */
        for (int j = 0; j < M; j++) {
            /* For all traversals of size (i,j) */
            /* window is 4 | constant */
            if (abs(i - j) > 0) break;
            for (int h = 0; h < TraversalsTable[i][j].size(); h++) {
                /* For all pairs of traversal h -- index 0 there are ids */
                termU.push_back(traversal_counter);
                termV.push_back(traversal_counter);
                Map_TraversalX_to_Curve.push_back(TraversalsTable[i][j][h][0][0][0]);
                Map_TraversalY_to_Curve.push_back(TraversalsTable[i][j][h][0][0][1]);
                for (int t = 1; t < TraversalsTable[i][j][h].size(); t++) {
                    vector<double> U = TraversalsTable[i][j][h][t][0];
                    vector<double> V = TraversalsTable[i][j][h][t][1];
                    for (int k = 0; k < G.size(); k++){
                        Unew = 0.0;
                        Vnew = 0.0;
                        for (int dim = 0; dim < d; dim++) {
                            Unew += G[k][dim] * U[dim];
                            Vnew += G[k][dim] * V[dim];
                        }
                        termU.push_back(Unew);
                        termV.push_back(Vnew);
                    }
                }
                Vectored_Traversals_X.push_back(termU);
                Vectored_Traversals_Y.push_back(termV);
                vector<double>().swap(termU);
                vector<double>().swap(termV);
                traversal_counter += 1;
            }
        }
        min_distance = new double [Vectored_Traversals_Y.size()];
        time = new double [Vectored_Traversals_Y.size()];
        nearest_neighbor = new int [Vectored_Traversals_Y.size()];
        /* init arrays */
        for (int i = 0; i < Vectored_Traversals_Y.size(); i++) {
            min_distance[i] = INT_MAX;
            nearest_neighbor[i] = -1;
            time[i] = 0;
        }
        cout << "Calling LSH ... " << Vectored_Traversals_X.size() << " " << Vectored_Traversals_Y.size() << endl;
        LSH(&Vectored_Traversals_X, &Vectored_Traversals_Y, k, L_vec, w, R, &R_neighbors, &min_distance, &time, &nearest_neighbor);
        /* index the nearest neighbor traversal to nearest neighbor id */
        vector<vector<int>> pair_ids;
        vector<int> ids;
        for (int n = 0; n < Vectored_Traversals_Y.size(); n++) {
            /* LSH returns the ids of the traversals, and we map them to the curve that they belong */
            ids.push_back((int)Map_TraversalY_to_Curve[Vectored_Traversals_Y[n][0]]);
            ids.push_back((int)Map_TraversalX_to_Curve[Vectored_Traversals_X[nearest_neighbor[n]][0]]);
            pair_ids.push_back(ids);
            traversal_neighbors.push_back(pair_ids);
            vector<int>().swap(ids);
            vector<vector<int>>().swap(pair_ids);
        }
        vector<vector<double>>().swap(Vectored_Traversals_X);
        vector<vector<double>>().swap(Vectored_Traversals_Y);
        vector<double>().swap(Map_TraversalX_to_Curve);
        vector<double>().swap(Map_TraversalY_to_Curve);
        delete[] min_distance;
        delete[] time;
        delete[] nearest_neighbor;
    }

    cout << "LSH traversals neighbors ready" << endl;

    double distance = 0.0;
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;
    double average_time = 0.0;

    min_distance = new double[searchset.size()];
    nearest_neighbor = new int[searchset.size()];
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = -1;
        nearest_neighbor[i] = -1;
    }

    /* for all lengths */
    for (int i = 0; i < traversal_neighbors.size(); i++){
        /* for all traversal curves */
        for (int j = 0; j < traversal_neighbors[i].size(); j++){
            cout <<traversal_neighbors[i][j][0]<< " | " << traversal_neighbors[i][j][1] << endl;
        }
    }
    /* ---- edit the results of the lsh ----*/
    double dtw_dist;
    /* for all lengths */
    for (int i = 0; i < traversal_neighbors.size(); i++){
        /* for all traversal curves */
        for (int j = 0; j < traversal_neighbors[i].size(); j++){
            if (traversal_neighbors[i][j][1] != -1 ) {
                dtw_dist = DTW(&dataset[traversal_neighbors[i][j][1]], &searchset[traversal_neighbors[i][j][0]]);
                if (min_distance[traversal_neighbors[i][j][1]] == -1) {
                    min_distance[traversal_neighbors[i][j][1]] = dtw_dist;
                    nearest_neighbor[traversal_neighbors[i][j][1]] = traversal_neighbors[i][j][0];
                } else if (distance < min_distance[traversal_neighbors[i][j][1]]) {
                    min_distance[traversal_neighbors[i][j][1]] = dtw_dist;
                    nearest_neighbor[traversal_neighbors[i][j][1]] = traversal_neighbors[i][j][0];
                }
             }
        }
    }
    /* ----- */

    /* Statistics available only when Brute Force is on */
    if (bfsearch == 'y' || bfsearch == 'Y') {
    /* compare with DTW the L different neighbors for every q*/
    /* Results for every curve query */
        int computations = 0;
        for (int q = 0; q < searchset.size(); q++) {
            curr_fraction = (double) min_distance[q] / TrueDistances[q];
            if (curr_fraction > max_af) max_af = curr_fraction;
            average_af += curr_fraction;
        }

        /* --- RESULTS --- */
        average_af = average_af / searchset.size();
        average_time = average_time / searchset.size();
        cout << "MAX Approximation Fraction (Grid/HyperCube Distance / True Distance) = " << max_af << endl;
        cout << "Average Approximation Fraction (Grid/HyperCube Distance / True Distance) = " << average_af << endl;
    }

    string Method = "LSH";
    ofstream neighbors_file;
    /* open file to dump all query results */
    neighbors_file.open("./output/curve_projection_lsh.txt");
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Query: " << i + 1 << endl;
        neighbors_file << "Method: LSH" << endl;
        neighbors_file << "HashFunction: " << Method << endl;
        if (nearest_neighbor[i] != -1) {
            neighbors_file << "Found Nearest Neighbor: " << nearest_neighbor[i] + 1 << endl;
        }else{
            neighbors_file << "Found Nearest Neighbor: Fail" << endl;
        }
        if (bfsearch == 'y' || bfsearch == 'Y') {
            neighbors_file << "True Nearest Neighbor: " << TrueNeighbors[i] + 1 << endl;
        }
        if (nearest_neighbor[i] != -1) {
            neighbors_file << "distanceFound: " << min_distance[i] << endl;
        }else{
            neighbors_file << "distanceFound: None" << endl;
        }
        if (bfsearch == 'y' || bfsearch == 'Y') {
            neighbors_file << "distanceTrue: " << TrueDistances[i] << endl;
        }
        neighbors_file << endl;
    }
    neighbors_file.close();

    cout << "Statistics are available only when Brute Force option is selected." << endl;

    /* clean remaining used memory */
    delete[] min_distance;
    delete[] nearest_neighbor;

    return 0;
}
