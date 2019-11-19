#include "Traversals.h"
#include "Helper_Functions.h"

void show_grid_lsh_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k_vec <int>         Number of hi function for construction of g function\n"
              << "\t-L_grid <int>         Number of Grids\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

void show_grid_bhc_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k_hypercube <int>         Dimension of Hypercube\n"
              << "\t-L_grid <int>         Number of Grids\n"
              << "\t-M <int>         Max allowed points to check\n"
              << "\t-probes <int>    Max allowed vertices to check\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

void shift_grid(vector<double>* orthogonal_grid, int delta, int d) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    /* t uniformly in [0,d) */
    uniform_real_distribution<double> distribution (0, d);
    for (int i = 0; i < d; i++) {
        orthogonal_grid->push_back(distribution(generator));
    }
}

void hash_curve(vector<vector<double>>* hashed_curve, vector<double*>* curve, vector<double>* orthogonal_grid, double delta, int d) {
    double* pi;
    vector<double> pi_new;
    vector<double> pi_old;
    /* first index has the id and the dimensions of the curve */
    pi = (*curve)[0];
    pi_new.push_back(pi[0]);
    pi_new.push_back(pi[1]);
    hashed_curve->push_back(pi_new);
    pi_old = pi_new;
    for (int i = 1; i < (*curve)[0][1]; i++) {
        pi = (*curve)[i];
        pi_new = arg_min(&pi, orthogonal_grid, delta, d);
        /* remove consecutive duplicates pi' from the hashed_curve */
        if ((pi_new[0] != pi_old[0]) || (pi_new[1] != pi_old[1])) {
            hashed_curve->push_back(pi_new);
        }
        pi_old = pi_new;
    }
}

void find_diagonal(vector<vector<int>>* v, int len1, int len2){
    double j = 0;
    double lamda = (double) len2 / len1;
    double stepi;
    if (len1>len2){
        stepi = (double) len2 / len1;
    } else{
        stepi = (double) len1 / len2;
    }
    vector<vector<int>> points;
    for (double i = 0; i < len1; i+=stepi) {
        vector<int> temp;
        int temp_i, temp_j;
        j = lamda * i;
        vector<int> pair1;
        pair1.push_back(floor(i));
        pair1.push_back(floor(j));
        points.push_back(pair1);
        temp_i = pair1[0] - 1;                               // element below
        temp_j = pair1[1];
        if((temp_i >= 0) && (temp_i < len1) && (temp_j >= 0) && (temp_j < len2)){
            temp.push_back(temp_i);
            temp.push_back(temp_j);
            points.push_back(temp);
            vector<int>().swap(temp);
        }
        temp_i = pair1[0] + 1;                               // element above
        temp_j = pair1[1];
        if((temp_i >= 0) && (temp_i < len1) && (temp_j >= 0) && (temp_j < len2)){
            temp.push_back(temp_i);
            temp.push_back(temp_j);
            points.push_back(temp);
            vector<int>().swap(temp);
        }
        if((ceil(j) != floor(j)) && (ceil(j) < len2)) {
            vector<int> pair2;
            pair2.push_back(floor(i));
            pair2.push_back(ceil(j));
            points.push_back(pair2);
            temp_i = pair2[0] - 1;                               // element below
            temp_j = pair2[1];
            if((temp_i >= 0) && (temp_i < len1) && (temp_j >= 0) && (temp_j < len2)){
                temp.push_back(temp_i);
                temp.push_back(temp_j);
                points.push_back(temp);
                vector<int>().swap(temp);
            }
            temp_i = pair2[0] + 1;                               // element above
            temp_j = pair2[1];
            if((temp_i >= 0) && (temp_i < len1) && (temp_j >= 0) && (temp_j < len2)){
                temp.push_back(temp_i);
                temp.push_back(temp_j);
                points.push_back(temp);
                vector<int>().swap(temp);
            }
        }
    }
    vector<int> temp_point;
    int found;
    for(int i = 0; i < points.size(); i++){
        found = 0;
        for ( int j = 0; j < v->size(); j++){
            if(points[i][0] == (*v)[j][0] && points[i][1] == (*v)[j][1]){
                found = 1;
                break;
            }
        }
        if ( found == 0 ){
            temp_point.push_back(points[i][0]);
            temp_point.push_back(points[i][1]);
            v->push_back(temp_point);
            vector<int>().swap(temp_point);
        }
    }
}

//points    : All points in diagonal +-1 cell
//i, j      : Current position of the robot (For the first call use 0,0)
//len1, len2: Dimensions of given the matrix
//traversal : The path traversed by robot till now (Vector to hold the path) */
void find_traversals(vector<vector<int>>* points, int i, int j, int len1, int len2, vector<vector<int>>* traversal, vector<vector<vector<int>>>* Traversals)
{
    vector<int> temp_step;
    // Reached the top of the matrix so we are left with
    // only option to move right


    // Add the current cell to the path being generated
    temp_step.push_back(i);
    temp_step.push_back(j);
    traversal->push_back(temp_step);
    vector<int>().swap(temp_step);

    if (i == len1 - 1)
    {
        for (int k = j+1; k < len2; k++){
            temp_step.push_back(i);
            temp_step.push_back(k);
            traversal->push_back(temp_step);
            vector<int>().swap(temp_step);
        }
        Traversals->push_back(*traversal);
        for (int k = j+1; k < len2; k++)
            traversal->erase(traversal->end());
        return;
    }

    // Reached the right corner of the matrix we are left with
    // only the upward movement.
    if (j == len2 - 1)
    {
        for (int k = i+1; k < len1; k++) {
            temp_step.push_back(k);
            temp_step.push_back(j);
            traversal->push_back(temp_step);
            vector<int>().swap(temp_step);
        }
        Traversals->push_back(*traversal);
        for (int k = i+1; k < len1; k++)
            traversal->erase(traversal->end());
        return;
    }

    // Print all the paths that are possible after moving up
    int found1 = 0;
    for (int k = 0; k < (*points).size(); k++){
        if (((*points)[k][0] == i+1) && ((*points)[k][1] == j)) {
            found1 = 1;
            break;
        }
    }
    if (found1 == 1) {
        find_traversals(points, i + 1, j, len1, len2, traversal,Traversals);
        traversal->erase(traversal->end());
    }

    // Print all the paths that are possible after moving right
    int found2 = 0;
    for (int k = 0; k < (*points).size(); k++){
        if (((*points)[k][0] == i) && ((*points)[k][1] == j+1)) {
            found2 = 1;
            break;
        }
    }
    if (found2 == 1) {
        find_traversals(points, i, j + 1, len1, len2, traversal,Traversals);
        traversal->erase(traversal->end());
    }

    //Print all the paths that are possible after moving diagonal
//    int found3 = 0;
//    for (int k = 0; k < (*points).size(); k++){
//        if (((*points)[k][0] == i+1) && ((*points)[k][1] == j+1)) {
//            found3 = 1;
//            break;
//        }
//    }
//    if (found3 == 1) {
//        find_traversals(points, i + 1, j + 1, len1, len2, traversal,Traversals);
//    }
}

void Relevant_Traversals(vector<vector<vector<int>>>* Traversals, int len1, int len2) {
    vector<vector<int>> points;

    find_diagonal(&points, len1, len2);
    vector<vector<int>> traversals;
    find_traversals(&points, 0, 0, len1, len2, &traversals, Traversals);
}
