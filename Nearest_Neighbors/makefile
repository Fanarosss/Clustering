all : lsh hypercube curves_grid_lsh curves_grid_hypercube curves_projection_lsh

curves_projection_lsh : HashTable.o Helper_Functions.o Traversals.o LSH_Functions.o LSH.o RandomProjections.o
	@echo "Compile projections curves ...";
	g++ -std=c++11 -I ./lib ./build/RandomProjections.o ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/Traversals.o ./build/HashTable.o -o ./build/curves_projection_lsh.x

RandomProjections.o :
	g++ -std=c++11 -I ./lib -c ./src/Projections_Curves/RandomProjections.cpp -o ./build/RandomProjections.o

curves_grid_hypercube : HashTable.o Helper_Functions.o Traversals.o BHC_Functions.o LSH_Functions.o LSH.o BHC.o Grid_Curves_bhc_main.o
	@echo "Compile grid curves hypercube ...";
	g++ -std=c++11 -I ./lib ./build/Grid_Curves_bhc_main.o ./build/BHC.o ./build/BHC_Functions.o ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/Traversals.o ./build/HashTable.o -o ./build/curves_grid_hypercube.x

curves_grid_lsh : HashTable.o Helper_Functions.o Traversals.o LSH_Functions.o LSH.o Grid_Curves_lsh_main.o
	@echo "Compile grid curves lsh ...";
	g++ -std=c++11 -I ./lib ./build/Grid_Curves_lsh_main.o ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/Traversals.o ./build/HashTable.o -o ./build/curves_grid_lsh.x

Grid_Curves_lsh_main.o :
	g++ -std=c++11 -I ./lib -c ./src/Grid_Curves/Curves.cpp -D _LSH_ -o ./build/Grid_Curves_lsh_main.o

Grid_Curves_bhc_main.o :
	g++ -std=c++11 -I ./lib -c ./src/Grid_Curves/Curves.cpp -D _BHC_ -o ./build/Grid_Curves_bhc_main.o

Traversals.o :
	g++ -std=c++11 -I ./lib -c ./src/Grid_Curves/Traversals.cpp -o ./build/Traversals.o

hypercube : HashTable.o Helper_Functions.o LSH_Functions.o BHC_Functions.o BHC.o BHC_main.o
	@echo "Compile hypercube ...";
	g++ -std=c++11 -I ./lib ./build/BHC_main.o ./build/BHC.o ./build/BHC_Functions.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/hypercube.x

BHC_main.o :
	g++ -std=c++11 -I ./lib -c ./src/Hypercube/BHC_main.cpp -o ./build/BHC_main.o

BHC.o :
	g++ -std=c++11 -I ./lib -c ./src/Hypercube/BHC.cpp -o ./build/BHC.o

BHC_Functions.o :
	g++ -std=c++11 -I ./lib -c ./src/Hypercube/BHC_Functions.cpp -o ./build/BHC_Functions.o

lsh : HashTable.o Helper_Functions.o LSH_Functions.o LSH.o LSH_main.o
	@echo "Compile lsh ...";
	g++ -std=c++11 -I ./lib ./build/LSH_main.o ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/lsh.x

LSH_main.o :
	g++ -std=c++11 -I ./lib -c ./src/LSH/LSH_main.cpp -o ./build/LSH_main.o

LSH.o :
	g++ -std=c++11 -I ./lib -c ./src/LSH/LSH.cpp -o ./build/LSH.o

LSH_Functions.o :
	g++ -std=c++11 -I ./lib -c ./src/LSH/LSH_Functions.cpp -o ./build/LSH_Functions.o

Helper_Functions.o :
	g++ -std=c++11 -I ./lib -c ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

HashTable.o :
	g++ -std=c++11 -I ./lib -c ./src/HashTable.cpp -o ./build/HashTable.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x
