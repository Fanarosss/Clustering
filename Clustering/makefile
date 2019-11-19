all : cluster

cluster : Helper_Functions.o Initializers.o Assigners.o Updaters.o Cluster.o Cluster_DataTypes.o main.o
	@echo "Compile Cluster ...";
	g++ -std=c++11 -I ./lib ./build/main.o ./build/Cluster_DataTypes.o ./build/Cluster.o ./build/Updaters.o ./build/Assigners.o ./build/Initializers.o ./build/Helper_Functions.o -o ./build/cluster.x

main.o :
	g++ -std=c++11 -I ./lib -c ./src/Clustering_main.cpp -o ./build/main.o

Cluster_DataTypes.o :
	g++ -std=c++11 -I ./lib -c ./src/Clusters/Cluster_DataTypes.cpp -o ./build/Cluster_DataTypes.o

Cluster.o :
	g++ -std=c++11 -I ./lib -c ./src/Clusters/Cluster.cpp -o ./build/Cluster.o

Updaters.o :
	g++ -std=c++11 -I ./lib -c ./src/Updaters/Updaters.cpp -o ./build/Updaters.o

Assigners.o :
	g++ -std=c++11 -I ./lib -c ./src/Assigners/Assigners.cpp -o ./build/Assigners.o

Initializers.o :
	g++ -std=c++11 -I ./lib -c ./src/Initializers/Initializers.cpp -o ./build/Initializers.o

Helper_Functions.o :
	g++ -std=c++11 -I ./lib -c ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x