#include <string>
#include <vector>
#include "Library.h"

using namespace std;

template <class Point>
class Updater {
protected:
    int K;
public:
    Updater(){}
    virtual void update(vector<vector<Point>>*, vector<int>**, vector<int>*) {}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class PAM : public Updater<Point> {
private:
    string name = "Partitioning Around Medoids (PAM)";
public:
    PAM(int K){this->K = K;}
    void update(vector<vector<Point>>*, vector<int>**, vector<int>*);
    string get_name();
};

template <class Point>
class MV_DTW : public Updater<Point> {
private:
    string name = "Mean Vector - DTW centroid Curve";
public:
    MV_DTW(int K){this->K = K;}
    void update(vector<vector<Point>>*, vector<int>**, vector<int>*);
    string get_name();
};