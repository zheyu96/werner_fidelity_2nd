#ifndef __GREEDY_H
#define __GREEDY_H

#include "../PathMethodBase/PathMethod.h"

using namespace std;

class Greedy : public PathMethod {
    vector<int> par;
    void BFS(int src, int dst);
    Path find_path(int src, int dst);
public:
    Greedy();
    ~Greedy();
    void build_paths(Graph _graph, vector<SDpair> _requests);
};

#endif