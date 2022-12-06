# parallel
Parallel algorithm with OpenMPI of Bayesian-based global search with Hooke–Jeeves local refinement for multi-objective optimization problems


Project should be compiled in build folder.
Use comand "cd build" to change directory to ./build.
Build result is in ./build/optimize. 

Use following commands to

build:
	
	make

clean up: 

	make clean

run (-np 2 specifies how many cores to use):
 
	mpirun -np 2 ./optimize