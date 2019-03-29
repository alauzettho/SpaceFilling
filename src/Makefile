CC=g++
OPT=-O3
OMP=-fopenmp
OPTC=-std=c++11
EXEC=spacefilling

all: $(EXEC)

spacefilling: main.o parser.o objectivefunction.o bfgs.o neldermead.o recuit.o
	$(CC) -o spacefilling main.o parser.o objectivefunction.o bfgs.o neldermead.o recuit.o $(OPT) $(OMP)

main.o: main.cpp
	$(CC) -o main.o -c main.cpp $(OPT) $(OMP)

parser.o: parser.cpp
	$(CC) -o parser.o -c parser.cpp $(OPT)

objectivefunction.o: objectivefunction.cpp
	$(CC) -o objectivefunction.o -c objectivefunction.cpp $(OPT) $(OPTC)

bfgs.o: bfgs.cpp
	$(CC) -o bfgs.o -c bfgs.cpp $(OPT)

neldermead.o: neldermead.cpp
	$(CC) -o neldermead.o -c neldermead.cpp $(OPT)

recuit.o: recuit.cpp
	$(CC) -o recuit.o -c recuit.cpp $(OPT) $(OMP)


clean:
	rm -rf *.o
	rm -rf $(EXEC)
	rm -rf ../output.txt

run:
	./spacefilling 2 16 RS MST 1