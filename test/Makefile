CPP=g++
CPPFLAGS=-O3 -Wall
LAPACK=-lblas -llapack

all: test_bsplines test_solver_banded test_interpolation test_interpolation_cubic

%:%.cpp *.h
	$(CC) $(CXXFLAGS) $*.cpp -o $* $(LAPACK)
