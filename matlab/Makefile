mex:
	mex -output rbfemmatrix rbfemmatrix.cpp CXXFLAGS='$$CXXFLAGS -g -Wall -pedantic -O3 '
oct:
	mkoctfile --mex -o rbfemmatrix rbfemmatrix.cpp
clean:
	rm rbfemmatrix.mex*

.DEFAULT_GOAL := mex

