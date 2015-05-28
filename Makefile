CC            = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/
LD_FLAGS        = -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
SAM             = samfire.o utilities_sam.o

samfire: $(SAM)
        $(CC) $(CC_FLAGS) $(SAM) -o samfire  $(LD_FLAGS)
samfire.o: samfire.cpp
        $(CC) $(CC_FLAGS) -c samfire.cpp
utilities_sam.o: utilities_sam.cpp
        $(CC) $(CC_FLAGS) -c utilities_sam.cpp
