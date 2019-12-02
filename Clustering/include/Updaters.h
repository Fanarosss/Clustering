#include <string>
#include <vector>
#include "Library.h"
#include "Database.h"

using namespace std;

template <class Point>
class Updater {
protected:
    int K;
public:
    Updater(){}
    virtual int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*) {return 0;}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class PAM : public Updater<Point> {
private:
    string name;
public:
    PAM(int);
    int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
    ~PAM();
};

template <class Point>
class MV_DTW : public Updater<Point> {
private:
    string name;
public:
    MV_DTW(int);
    int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
    ~MV_DTW();
};