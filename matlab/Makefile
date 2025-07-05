mex:
	mex -output rbfemmatrix rbfemmatrix.cpp CXXFLAGS='$$CXXFLAGS -g -Wall -pedantic -O3 '
oct:
	mkoctfile --mex -o rbfemmatrix rbfemmatrix.cpp
clean:
	rm rbfemmatrix.mex*

# derived the astyle settings from https://github.com/nlohmann/json
pretty:
	astyle \
	    --style=attach \
	    --indent=spaces=4 \
	    --indent-modifiers \
	    --indent-switches \
	    --indent-preproc-block \
	    --indent-preproc-define \
	    --indent-col1-comments \
	    --pad-oper \
	    --pad-header \
	    --align-pointer=type \
	    --align-reference=type \
	    --add-brackets \
	    --convert-tabs \
	    --close-templates \
	    --lineend=linux \
	    --preserve-date \
	    --suffix=none \
	    --formatted \
	    --break-blocks \
	   "*.h" "*.cpp"

.DEFAULT_GOAL := mex

