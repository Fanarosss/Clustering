all : cluster

cluster : Helper_Functions.o DistanceDatabase.o Initializers.o Assigners.o Updaters.o HashTable.o LSH_DataTypes.o LSH_Functions.o LSH.o Cluster.o Cluster_DataTypes.o main.o
	@echo "Compile Cluster ...";
	g++ -std=c++11 -I ./include ./build/main.o ./build/Cluster_DataTypes.o ./build/Cluster.o ./build/LSH.o ./build/LSH_Functions.o ./build/LSH_DataTypes.o ./build/HashTable.o ./build/Updaters.o ./build/Assigners.o ./build/Initializers.o ./build/DistanceDatabase.o ./build/Helper_Functions.o -o ./build/cluster.x

main.o :
	g++ -std=c++11 -I ./include -c ./src/Clustering_main.cpp -o ./build/main.o

Cluster_DataTypes.o :
	g++ -std=c++11 -I ./include -c ./src/Clusters/Cluster_DataTypes.cpp -o ./build/Cluster_DataTypes.o

Cluster.o :
	g++ -std=c++11 -I ./include -c ./src/Clusters/Cluster.cpp -o ./build/Cluster.o

LSH.o :
	g++ -std=c++11 -I ./include -c ./src/LSH/LSH.cpp -o ./build/LSH.o

LSH_Functions.o :
	g++ -std=c++11 -I ./include -c ./src/LSH/LSH_Functions.cpp -o ./build/LSH_Functions.o

LSH_DataTypes.o :
	g++ -std=c++11 -I ./include -c ./src/LSH/LSH_DataTypes.cpp -o ./build/LSH_DataTypes.o

HashTable.o :
	g++ -std=c++11 -I ./include -c ./src/LSH/HashTable.cpp -o ./build/HashTable.o

Updaters.o :
	g++ -std=c++11 -I ./include -c ./src/Updaters/Updaters.cpp -o ./build/Updaters.o

Assigners.o :
	g++ -std=c++11 -I ./include -c ./src/Assigners/Assigners.cpp -o ./build/Assigners.o

Initializers.o :
	g++ -std=c++11 -I ./include -c ./src/Initializers/Initializers.cpp -o ./build/Initializers.o

DistanceDatabase.o :
	g++ -std=c++11 -I ./include -c ./src/Database.cpp -o ./build/DistanceDatabase.o

Helper_Functions.o :
	g++ -std=c++11 -I ./include -c ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x