IDIR=include
GCC=g++
LIBS=-fopenmp -lmpfr -lgmp 
CPPFLAGS= -I$(IDIR)   -Wall -fopenmp -std=c++17 -O3 -g -rdynamic -march=native -D EIGEN_DONT_PARALLELIZE -D EIGEN_DONT_VECTORIZE
SOURCES=$(wildcard src/*.cpp)
OBJ_1=$(patsubst %.cpp,%.o,$(SOURCES))
OBJ=$(patsubst src/%, obj/%, $(OBJ_1))
DEP1= $(patsubst %.cpp,%.d,$(SOURCES))
DEP=$(patsubst src/%, dep/%, $(DEP1))
$(shell mkdir -p obj)
dep/%.d: src/%.cpp
	@set -e; rm -f $@; \
	$(GCC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

obj/%.o: src/%.cpp
	$(GCC) -c  $(CPPFLAGS) -o $@ $<


HLT_spallate: $(OBJ)
	$(GCC)  -o $@ $^ $(CPPFLAGS) -L $(LIB_PY) $(LIBS)



.PHONY: clean

clean :
	rm -f HLT_spallate $(OBJ) $(DEP)
