CC = g++
CFLAGS = -Wall -O3 -std=c++11
LFLAGS = 
OBJECTS = obj/main.o obj/correlationMatrix.o obj/dataReader.o obj/heatMap.o obj/k_means.o obj/stats.o obj/unsupervisedNormalization.o obj/utilities.o obj/lodepng.o
EXECUTABLE = CancerClustering

all : $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $(EXECUTABLE)

obj/main.o : src/main.cpp
	$(CC) $(CFLAGS) -c src/main.cpp -o obj/main.o

obj/correlationMatrix.o : src/correlationMatrix.cpp
	$(CC) $(CFLAGS) -c src/correlationMatrix.cpp -o obj/correlationMatrix.o

obj/dataReader.o : src/dataReader.cpp
	$(CC) $(CFLAGS) -c src/dataReader.cpp -o obj/dataReader.o
	
obj/heatMap.o : src/heatMap.cpp
	$(CC) $(CFLAGS) -c src/heatMap.cpp -o obj/heatMap.o

obj/k_means.o : src/k_means.cpp
	$(CC) $(CFLAGS) -c src/k_means.cpp -o obj/k_means.o
	
obj/stats.o : src/stats.cpp
	$(CC) $(CFLAGS) -c src/stats.cpp -o obj/stats.o
	
obj/unsupervisedNormalization.o : src/unsupervisedNormalization.cpp
	$(CC) $(CFLAGS) -c src/unsupervisedNormalization.cpp -o obj/unsupervisedNormalization.o
	
obj/utilities.o : src/utilities.cpp
	$(CC) $(CFLAGS) -c src/utilities.cpp -o obj/utilities.o
	
obj/lodepng.o : src/lodePNG/lodepng.cpp
	$(CC) $(CFLAGS) -c src/lodePNG/lodepng.cpp -o obj/lodepng.o
	
clean:
	rm *o hello