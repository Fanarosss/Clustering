GOOGLE_TEST_LIB = gtest
GOOGLE_TEST_INCLUDE = /usr/local/include

G++ = g++
G++_GTEST_FLAGS = -c -Wall -I $(GOOGLE_TEST_INCLUDE) -I ./include
LD_GTEST_FLAGS = -L /usr/local/lib -l $(GOOGLE_TEST_LIB) -l pthread
G++_FLAGS = -c -ggdb3 -I ./include

all : cluster

tests: gtest_main.o gtest_modulo.o gtest_database.o DistanceDatabase.o Helper_Functions.o
	g++ -std=c++11 -o ./build/gtest_modulo.o ./build/gtest_database.o ./build/DistanceDatabase.o ./build/Helper_Functions.o ./build/gtest_main.o $(LD_GTEST_FLAGS) -o ./build/execute_tests.x

gtest_main.o :
	g++ -std=c++11 $(G++_GTEST_FLAGS) ./tests/gtest_main.cpp -o ./build/gtest_main.o

gtest_modulo.o :
	g++ -std=c++11 $(G++_GTEST_FLAGS) ./tests/gtest_modulo.cpp -o ./build/gtest_modulo.o

gtest_database.o :
	g++ -std=c++11 $(G++_GTEST_FLAGS) ./tests/gtest_database.cpp -o ./build/gtest_database.o

cluster : Helper_Functions.o DistanceDatabase.o Initializers.o Assigners.o Updaters.o HashTable.o LSH_DataTypes.o LSH_Functions.o LSH.o Cluster.o Cluster_DataTypes.o main.o
	@echo "Compile Cluster ...";
	g++ -std=c++11 ./build/main.o ./build/Cluster_DataTypes.o ./build/Cluster.o ./build/LSH.o ./build/LSH_Functions.o ./build/LSH_DataTypes.o ./build/HashTable.o ./build/Updaters.o ./build/Assigners.o ./build/Initializers.o ./build/DistanceDatabase.o ./build/Helper_Functions.o -o ./build/cluster.x

main.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Clustering_main.cpp -o ./build/main.o

Cluster_DataTypes.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Clusters/Cluster_DataTypes.cpp -o ./build/Cluster_DataTypes.o

Cluster.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Clusters/Cluster.cpp -o ./build/Cluster.o

LSH.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/LSH/LSH.cpp -o ./build/LSH.o

LSH_Functions.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/LSH/LSH_Functions.cpp -o ./build/LSH_Functions.o

LSH_DataTypes.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/LSH/LSH_DataTypes.cpp -o ./build/LSH_DataTypes.o

HashTable.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/LSH/HashTable.cpp -o ./build/HashTable.o

Updaters.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Updaters/Updaters.cpp -o ./build/Updaters.o

Assigners.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Assigners/Assigners.cpp -o ./build/Assigners.o

Initializers.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Initializers/Initializers.cpp -o ./build/Initializers.o

DistanceDatabase.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Database.cpp -o ./build/DistanceDatabase.o

Helper_Functions.o :
	g++ -std=c++11 $(G++_FLAGS) ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x