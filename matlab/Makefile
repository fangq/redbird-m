mex:
	mex -output rbfemmatrix femdiffusion.cpp CXXFLAGS='$$CXXFLAGS -g'
oct:
	mkoctfile --mex -o rbfemmatrix femdiffusion.cpp
clean:
	rm rbfemmatrix.mex*

.DEFAULT_GOAL := mex

