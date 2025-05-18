
dependency = Network/Node/Node.o Network/Shape/Shape.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o config.o
algorithm = Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyAlgo1/MyAlgo1.o Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo4/MyAlgo4.o Algorithm/MyAlgo5/MyAlgo5.o Algorithm/MyAlgo6/MyAlgo6.o
algorithm2 = Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo1_time/MyAlgo1_time.o
routing = Network/PathMethod/PathMethodBase/PathMethod.o Network/PathMethod/Greedy/Greedy.o Network/PathMethod/QCAST/QCAST.o Network/PathMethod/REPS/REPS.o
gcc_flag = -g -Wall -Wextra -Os -Ofast --std=c++17

gurobi_path = /home/oba/danger/gurobi951/linux64
gurobi_flag = -I/home/oba/danger/gurobi951/linux64/include/ -L/home/oba/danger/gurobi951/linux64/lib/ -lgurobi_c++ -lgurobi95

all: a.out

runtime: runtime.cpp $(dependency) $(algorithm2) $(routing)
	g++ -fopenmp $(gcc_flag) runtime.cpp $(dependency) $(algorithm2) -o runtime

a.out: main.cpp $(dependency) $(algorithm) $(routing)
	g++ -fopenmp $(gcc_flag) main.cpp $(dependency) $(algorithm) $(routing) -I/home/oba/danger/gurobi951/linux64/include/ -L/home/oba/danger/gurobi951/linux64/lib/ -lgurobi_c++ -lgurobi95

cdf: cdf.cpp $(dependency) $(algorithm) $(routing)
	g++ -fopenmp $(gcc_flag) cdf.cpp $(dependency) $(algorithm) -o cdf

Network/Node/Node.o: Network/Node/Node.cpp Network/Node/Node.h config.h
	g++ -c $(gcc_flag) Network/Node/Node.cpp Network/Node/Node.h config.h
	mv Node.o Network/Node/Node.o

Network/Shape/Shape.o: Network/Shape/Shape.cpp Network/Shape/Shape.h config.h
	g++ -c $(gcc_flag) Network/Shape/Shape.cpp Network/Shape/Shape.h config.h
	mv Shape.o Network/Shape/Shape.o

Network/Graph/Graph.o: Network/Graph/Graph.cpp Network/Graph/Graph.h config.h
	g++ -c $(gcc_flag) Network/Graph/Graph.cpp Network/Graph/Graph.h config.h
	mv Graph.o Network/Graph/Graph.o

Algorithm/AlgorithmBase/AlgorithmBase.o: Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	g++ -c $(gcc_flag) Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	mv AlgorithmBase.o Algorithm/AlgorithmBase/AlgorithmBase.o

Algorithm/MyAlgo1/MyAlgo1.o: Algorithm/MyAlgo1/MyAlgo1.h Algorithm/MyAlgo1/MyAlgo1.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo1/MyAlgo1.h Algorithm/MyAlgo1/MyAlgo1.cpp config.h
	mv MyAlgo1.o Algorithm/MyAlgo1/MyAlgo1.o

Algorithm/MyAlgo1_time/MyAlgo1_time.o: Algorithm/MyAlgo1_time/MyAlgo1_time.h Algorithm/MyAlgo1_time/MyAlgo1_time.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo1_time/MyAlgo1_time.h Algorithm/MyAlgo1_time/MyAlgo1_time.cpp config.h
	mv MyAlgo1_time.o Algorithm/MyAlgo1_time/MyAlgo1_time.o

Algorithm/MyAlgo2/MyAlgo2.o: Algorithm/MyAlgo2/MyAlgo2.h Algorithm/MyAlgo2/MyAlgo2.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo2/MyAlgo2.h Algorithm/MyAlgo2/MyAlgo2.cpp config.h
	mv MyAlgo2.o Algorithm/MyAlgo2/MyAlgo2.o

Algorithm/MyAlgo3/MyAlgo3.o: Algorithm/MyAlgo3/MyAlgo3.h Algorithm/MyAlgo3/MyAlgo3.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo3/MyAlgo3.h Algorithm/MyAlgo3/MyAlgo3.cpp config.h
	mv MyAlgo3.o Algorithm/MyAlgo3/MyAlgo3.o

Algorithm/MyAlgo4/MyAlgo4.o: Algorithm/MyAlgo4/MyAlgo4.h Algorithm/MyAlgo4/MyAlgo4.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo4/MyAlgo4.h Algorithm/MyAlgo4/MyAlgo4.cpp config.h
	mv MyAlgo4.o Algorithm/MyAlgo4/MyAlgo4.o

Algorithm/MyAlgo5/MyAlgo5.o: Algorithm/MyAlgo5/MyAlgo5.h Algorithm/MyAlgo5/MyAlgo5.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo5/MyAlgo5.h Algorithm/MyAlgo5/MyAlgo5.cpp config.h
	mv MyAlgo5.o Algorithm/MyAlgo5/MyAlgo5.o

Algorithm/MyAlgo6/MyAlgo6.o: Algorithm/MyAlgo6/MyAlgo6.h Algorithm/MyAlgo6/MyAlgo6.cpp config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo6/MyAlgo6.h Algorithm/MyAlgo6/MyAlgo6.cpp config.h
	mv MyAlgo6.o Algorithm/MyAlgo6/MyAlgo6.o

Network/PathMethod/Greedy/Greedy.o: Network/PathMethod/Greedy/Greedy.h Network/PathMethod/Greedy/Greedy.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/Greedy/Greedy.h Network/PathMethod/Greedy/Greedy.cpp config.h
	mv Greedy.o Network/PathMethod/Greedy/Greedy.o

Network/PathMethod/QCAST/QCAST.o: Network/PathMethod/QCAST/QCAST.h Network/PathMethod/QCAST/QCAST.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/QCAST/QCAST.h Network/PathMethod/QCAST/QCAST.cpp config.h
	mv QCAST.o Network/PathMethod/QCAST/QCAST.o

Network/PathMethod/REPS/REPS.o: Network/PathMethod/REPS/REPS.h Network/PathMethod/REPS/REPS.cpp config.h
	g++ -c --static $(gcc_flag) Network/PathMethod/REPS/REPS.h Network/PathMethod/REPS/REPS.cpp config.h -I/home/oba/danger/gurobi951/linux64/include/ -L/home/oba/danger/gurobi951/linux64/lib/ -lgurobi_c++ -lgurobi95
	mv REPS.o Network/PathMethod/REPS/REPS.o

Network/PathMethod/PathMethodBase/PathMethod.o: Network/PathMethod/PathMethodBase/PathMethod.h Network/PathMethod/PathMethodBase/PathMethod.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/PathMethodBase/PathMethod.h Network/PathMethod/PathMethodBase/PathMethod.cpp config.h
	mv PathMethod.o Network/PathMethod/PathMethodBase/PathMethod.o

config.o: config.h config.cpp
	g++ -c $(gcc_flag) config.h config.cpp

clean:
	rm -f a.out
	find . -type f -name '*.o' -delete
	find . -type f -name '*.h.gch' -delete
	rm -f ../data/ans/*.ans
	rm -f ../data/log/*
	rm -f ../data/input/*
	rm cdf
	rm runtime