#include <string>
#include <vector>

using namespace std;

template <class Point>
class Initializer {
protected:
    int K;
public:
    Initializer(){}
    virtual vector<int>* init(vector<vector<Point>>*) {return NULL;}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class Random_Selection : public Initializer<Point> {
private:
    string name = "Random Selection";
public:
    Random_Selection(int K){this->K = K;}
    vector<int>* init(vector<vector<Point>>*);
    string get_name();
};

template <class Point>
class KMeans_plusplus : public Initializer<Point> {
private:
    string name = "K-Means++";
public:
    KMeans_plusplus(int K){this->K = K;}
    vector<int>* init(vector<vector<Point>>*);
    string get_name();
};