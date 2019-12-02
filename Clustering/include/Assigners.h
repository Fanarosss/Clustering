#include <string>
#include <vector>
#include "Library.h"

using namespace std;

template <class Point>
class Assigner {
protected:
    int K;
    int Grids;
    int L;
    int k;
public:
    Assigner(){}
    virtual vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>) {return NULL;}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class Lloyd_assignment : public Assigner<Point> {
private:
    string name;
public:
    Lloyd_assignment(int, int, int, int);
    vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>);
    string get_name();
    ~Lloyd_assignment();
};

template <class Point>
class Inverse_assignment : public Assigner<Point> {
private:
    string name;
public:
    Inverse_assignment(int, int, int, int);
    vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>);
    string get_name();
    ~Inverse_assignment();
};