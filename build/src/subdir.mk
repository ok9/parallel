OBJS += \
./src/Solution.o \
./src/Function.o \
./src/OptimizedFunction.o \
./src/SurogateFunction.o \
./src/SurogateOneCriteriaFunction.o \
./src/ParallelHybridAlgorithm.o \
./src/HookeJeeves.o \
./src/Hypercube.o \
./src/AllTestFunctions.o \
./src/TestParallelHybridAlgorithm.o





src/%.o: ../src/%.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -o "$@" -c "$<"

