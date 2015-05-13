CC = g++
CFLAGS = -Wall -O3 -std=c++11
LFLAGS = 

ODIR = obj
SDIR = src
TDIR = tests

_OBJ_CORE = main.o correlationMatrix.o dataReader.o heatMap.o k_means.o stats.o unsupervisedNormalization.o utilities.o lodePNG/lodepng.o
OBJ_CORE = $(patsubst %,$(ODIR)/%,$(_OBJ_CORE))

_OBJ_TEST = correlationMatrix_test.o dataReader_test.o general_test.o k_means_test.o stats_test.o unsupervisedNormalization_test.o utilities_test.o lodePNG_test.o
OBJ_TEST = $(patsubst %,$(ODIR)/$(TDIR)/%,$(_OBJ_TEST))

_HEADERS_CORE = correlationMatrix.hpp dataReader.hpp heatMap.hpp k_means.hpp stats.hpp unsupervisedNormalization.hpp utilities.hpp lodePNG/lodepng.h
HEADERS_CORE = $(patsubst %,$(SDIR)/%,$(_HEADERS_CORE))

EXECUTABLE = CancerClustering

all: ${EXECUTABLE}

${EXECUTABLE}: ${OBJ_CORE} ${OBJ_TEST} ${HEADERS_CORE}
	$(CC) $(LFLAGS) ${OBJ_CORE} ${OBJ_TEST} -o $(EXECUTABLE)

$(ODIR)/%.o: $(SDIR)/%.cpp ${HEADERS_CORE}
	$(CC) $(CFLAGS) -c $< -o $@	
	
$(ODIR)/$(TDIR)/%.o: $(SDIR)/$(TDIR)/%.cpp $(SDIR)/$(TDIR)/%.hpp ${HEADERS_CORE}
	$(CC) $(CFLAGS) -c $< -o $@	

.PHONY : clean

clean:
	rm ${EXECUTABLE} ${OBJ_CORE} ${OBJ_TEST}