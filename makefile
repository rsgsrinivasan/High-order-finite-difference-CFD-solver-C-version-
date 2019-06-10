SIM: Main.o characteristic_plane.o Flow_Initialisation.o inputfile_preprocessor.o Grid_Transformation.o Neighbour_search.o node_locations.o presolver.o Roe_averaging.o Solver.o Subdomain_neighbours.o Tecplot_output_files.o viscous_terms.o
	mpic++ -O3   -g  -pedantic -lm -std=c++0x   Main.o characteristic_plane.o Flow_Initialisation.o inputfile_preprocessor.o Grid_Transformation.o Neighbour_search.o node_locations.o presolver.o Roe_averaging.o Solver.o Subdomain_neighbours.o Tecplot_output_files.o viscous_terms.o -o SIM
Main.o: Main.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Main.cpp
characteristic_plane.o: characteristic_plane.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    characteristic_plane.cpp
Flow_Initialisation.o: Flow_Initialisation.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Flow_Initialisation.cpp	
Grid_Transformation.o: Grid_Transformation.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Grid_Transformation.cpp
Neighbour_search.o: Neighbour_search.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Neighbour_search.cpp
node_locations.o: node_locations.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    node_locations.cpp
presolver.o: presolver.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    presolver.cpp
Roe_averaging.o: Roe_averaging.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Roe_averaging.cpp
Solver.o: Solver.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Solver.cpp
Subdomain_neighbours.o: Subdomain_neighbours.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Subdomain_neighbours.cpp
Tecplot_output_files.o: Tecplot_output_files.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    Tecplot_output_files.cpp
viscous_terms.o: viscous_terms.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    viscous_terms.cpp
inputfile_preprocessor.o: inputfile_preprocessor.cpp functions.h
	mpic++ -O3   -g  -pedantic -lm -std=c++0x     -c    inputfile_preprocessor.cpp



clean:
		rm -rf *.o SIM gnu_plot.txt *.txt metis*
cleanall:
		rm -rf *.o SIM *.plt gnu_plot.txt *.txt metis* restart_file.neu node_neighbour.neu elem_neighbour.neu restart_file_1.neu restart_file_2.neu *.dat
cleanr:
		rm -rf *.o SIM *.plt gnu_plot.txt *.txt metis* restart_file.neu
