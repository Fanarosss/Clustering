#include "Library.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "Traversals.h"

using namespace std;

void Grid_Vectorization(double delta, int d, vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, vector<vector<double>>* data_vectored_curves, vector<vector<double>>* search_vectored_curves) {

    /* variable decl */
    double max_element = 0.0;
    int max_points = 0, elements = 0;
    /* orthogonal grid of size d */
    vector<double> orthogonal_grid;
    /* vector for hashed curves */
    vector<vector<vector<double>>> hashed_curves;
    vector<vector<double>> temp_hash;
    /* temp */
    vector<double> curve;

    /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
    shift_grid(&orthogonal_grid, delta, d);

    /* ------------------ DATA SET hashing ----------------- */
    /* hash all dataset curves */
    for (int i = 0; i < dataset->size(); i++) {
        hash_curve(&temp_hash, &(*dataset)[i], &orthogonal_grid, delta, d);
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
        data_vectored_curves->push_back(curve);
        vector<double>().swap(curve);
    }

    /* ------------------ SEARCH SET hashing ----------------- */
    /* clean hashed curves */
    vector<vector<vector<double>>>().swap(hashed_curves);
    /* end of cleaning */

    /* hash all curves */
    for (int i = 0; i < searchset->size(); i++) {
        hash_curve(&temp_hash, &(*searchset)[i], &orthogonal_grid, delta, d);
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
        search_vectored_curves->push_back(curve);
        vector<double>().swap(curve);
    }

    /* clear no longer used vectors for memory optimizations */
    vector<double>().swap(orthogonal_grid);
    /* clean hashed curves */
    vector<vector<vector<double>>>().swap(hashed_curves);
    /* end of cleaning */

    /* ----------------- PADDING for both sets ------------ */
    /* pad special number > max coord */
    for (int i = 0; i < data_vectored_curves->size(); i++) {
        while ((*data_vectored_curves)[i].size() < max_points) {
            (*data_vectored_curves)[i].push_back(2 * max_element);
        }
    }
    for (int i = 0; i < search_vectored_curves->size(); i++) {
        while ((*search_vectored_curves)[i].size() < max_points) {
            (*search_vectored_curves)[i].push_back(2 * max_element);
        }
    }
    /* ----------- end of padding ------------ */
    return;
}

