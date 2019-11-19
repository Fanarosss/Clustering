#include "HashTable.h"
#include "Helper_Functions.h"

using namespace std;

template <class data>
HashTable<data>::HashTable(int TableSize) {
    pTable = new vector<vector<data>> [TableSize];
    Size = TableSize;
}

template <class data>
int HashTable<data>::Hash(int HashCode) {
    return HashCode % Size;
}

template <class data>
void HashTable<data>::Insert(int HashCode, vector<data> point){
    int hashed = Hash(HashCode);
    pTable[hashed].push_back(point);
}

template <class data>
vector<vector<data>>* HashTable<data>::Search_Neighbors(int HashCode){
    int hashed = Hash(HashCode);
    return &pTable[hashed];
}

template <class data>
HashTable<data>::~HashTable(){
    delete[] pTable;
    return;
}

template class HashTable<int>;
template class HashTable<double>;
