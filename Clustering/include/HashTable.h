#include "Library.h"

using namespace std;

template <class data>
class HashTable {
private:
    vector<vector<data>>* pTable;
    int Size;
public:
    HashTable(int);
    int Hash(int);
    void Insert(int, vector<data>);
    vector<vector<data>>* Search_Neighbors(int);
    ~HashTable();
};
