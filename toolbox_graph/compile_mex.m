% compile mex file for dijkstra

mex mex/perform_dijkstra_propagation.cpp mex/fheap/fib.cpp 
mex mex/dijkstra.cpp -o perform_dijkstra_fast 