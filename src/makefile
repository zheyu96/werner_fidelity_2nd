# --- 路徑與變數設定 ---
dependency = Network/Node/Node.o Network/Shape/Shape.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o config.o
algorithm = Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyAlgo1/MyAlgo1.o Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo4/MyAlgo4.o \
            Algorithm/MyAlgo5/MyAlgo5.o Algorithm/MyAlgo6/MyAlgo6.o Algorithm/WernerAlgo/WernerAlgo.o \
            Algorithm/WernerAlgo2/WernerAlgo2.o Algorithm/WernerAlgo_UB/WernerAlgo_UB.o
algorithm2 = Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo1_time/MyAlgo1_time.o
algorithm_zfa = Algorithm/WernerAlgo_time/WernerAlgo_time.o 
routing = Network/PathMethod/PathMethodBase/PathMethod.o Network/PathMethod/Greedy/Greedy.o \
          Network/PathMethod/QCAST/QCAST.o Network/PathMethod/REPS/REPS.o

# 編譯旗標
gcc_flag = -g -Wall -Wextra -Os -Ofast --std=c++17

# Gurobi 設定
gurobi_path = /home/oba/danger/gurobi951/linux64
gurobi_flag = -I$(gurobi_path)/include/ -L$(gurobi_path)/lib/ -lgurobi_c++ -lgurobi95

# --- 主要執行目標 ---

all: a.out

# 專門量測 WernerAlgo_time 的目標
zfa_run: ZFA_runtime.cpp $(dependency) $(algorithm_zfa) $(routing)
	g++ -fopenmp $(gcc_flag) ZFA_runtime.cpp $(dependency) $(algorithm_zfa) $(routing) $(gurobi_flag) -o zfa_run
	./zfa_run

# 舊有的 runtime 目標 (修正：加上 $(routing) 與 $(gurobi_flag) 以免報錯)
runtime: runtime.cpp $(dependency) $(algorithm2) $(routing)
	g++ -fopenmp $(gcc_flag) runtime.cpp $(dependency) $(algorithm2) $(routing) $(gurobi_flag) -o runtime

# 預設的主程式目標
a.out: main.cpp $(dependency) $(algorithm) $(routing)
	g++ -fopenmp $(gcc_flag) main.cpp $(dependency) $(algorithm) $(routing) $(gurobi_flag) -o a.out
	./a.out

# --- 各個物件檔案 (Object Files) 編譯規則 ---

Network/Node/Node.o: Network/Node/Node.cpp Network/Node/Node.h config.h
	g++ -c $(gcc_flag) Network/Node/Node.cpp -o Network/Node/Node.o

Network/Shape/Shape.o: Network/Shape/Shape.cpp Network/Shape/Shape.h config.h
	g++ -c $(gcc_flag) Network/Shape/Shape.cpp -o Network/Shape/Shape.o

Network/Graph/Graph.o: Network/Graph/Graph.cpp Network/Graph/Graph.h config.h
	g++ -c $(gcc_flag) Network/Graph/Graph.cpp -o Network/Graph/Graph.o

Algorithm/AlgorithmBase/AlgorithmBase.o: Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	g++ -c $(gcc_flag) Algorithm/AlgorithmBase/AlgorithmBase.cpp -o Algorithm/AlgorithmBase/AlgorithmBase.o

# WernerAlgo_time 規則
Algorithm/WernerAlgo_time/WernerAlgo_time.o: Algorithm/WernerAlgo_time/WernerAlgo_time.h Algorithm/WernerAlgo_time/WernerAlgo_time.cpp config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo_time/WernerAlgo_time.cpp -o Algorithm/WernerAlgo_time/WernerAlgo_time.o

# 路由相關規則 (注意 REPS 需要連結 Gurobi)
Network/PathMethod/REPS/REPS.o: Network/PathMethod/REPS/REPS.h Network/PathMethod/REPS/REPS.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/REPS/REPS.cpp $(gurobi_flag) -o Network/PathMethod/REPS/REPS.o

Network/PathMethod/Greedy/Greedy.o: Network/PathMethod/Greedy/Greedy.h Network/PathMethod/Greedy/Greedy.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/Greedy/Greedy.cpp -o Network/PathMethod/Greedy/Greedy.o

Network/PathMethod/QCAST/QCAST.o: Network/PathMethod/QCAST/QCAST.h Network/PathMethod/QCAST/QCAST.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/QCAST/QCAST.cpp -o Network/PathMethod/QCAST/QCAST.o

Network/PathMethod/PathMethodBase/PathMethod.o: Network/PathMethod/PathMethodBase/PathMethod.h Network/PathMethod/PathMethodBase/PathMethod.cpp config.h
	g++ -c $(gcc_flag) Network/PathMethod/PathMethodBase/PathMethod.cpp -o Network/PathMethod/PathMethodBase/PathMethod.o

# 其餘演算法編譯規則 (範例：MyAlgo)
Algorithm/MyAlgo1/MyAlgo1.o: Algorithm/MyAlgo1/MyAlgo1.cpp Algorithm/MyAlgo1/MyAlgo1.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo1/MyAlgo1.cpp -o Algorithm/MyAlgo1/MyAlgo1.o

config.o: config.h config.cpp
	g++ -c $(gcc_flag) config.cpp -o config.o

# --- 清除工具 ---
clean:
	rm -f a.out zfa_run runtime
	find . -type f -name '*.o' -delete
	find . -type f -name '*.h.gch' -delete
	rm -f ../data/ans/*.ans
	rm -f ../data/log/*
	rm -f ../data/input/*

.PHONY: clean all zfa_run runtime