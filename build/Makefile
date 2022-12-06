FLAGS := -O2 -g -Wall
CXX   :=  mpicxx

EXECUTABLE := optimize
ARGS := 

-include src/subdir.mk

all: surrogates

surrogates: $(OBJS)
	mpicxx  -o $(EXECUTABLE) $(OBJS) $(LIBPATH) $(LIBS) -static-libstdc++ -static-libgcc
clean: 
	rm $(OBJS) $(EXECUTABLE)
	



