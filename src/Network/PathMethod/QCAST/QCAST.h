#ifndef __QCAST_H
#define __QCAST_H

#include "../PathMethodBase/PathMethod.h"
using namespace std;

/* 
找 major path 的方法: Greedy-EDA(G-EDA)
Step1)  用 Extended Dijkstra's algorithm 根據權重 EXT 找出每個 S-D pair 的最佳 path
Step2)  選最大 EXT 的 path，分配資源後更新為 residual graph
重複 step 1, 2
直到無法再找到 path，或找到的 path 超過一定數量
設定 path length 的最大值，超過就直接 pop 掉
(先不用recovery path)
*/ 
class QCAST : public PathMethod {
    double EXT(Path path, int w);
    map<pair<int, int>, double> combination;
    double C(int n, int m);
    pair<double, Path> find_best_path(int src, int dst);
public:
    QCAST();
    ~QCAST();
    void build_paths(Graph graph, vector<SDpair> requests);
};

#endif