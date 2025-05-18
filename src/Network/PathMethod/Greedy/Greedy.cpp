#include "Greedy.h"

Greedy::Greedy() {
    method_name = "Greedy";
}
Greedy::~Greedy() {}

void Greedy::build_paths(Graph _graph, vector<SDpair> _requests) {
    graph = _graph;
    requests = _requests;
    paths.clear();
    while(1) {
        vector<Path> candidate_path;
        for(auto [src, dst] : requests) {
            Path path = find_path(src, dst);
            if(path.empty()) continue;
            candidate_path.push_back(path);
        }

        if(candidate_path.empty()) break;

        sort(candidate_path.begin(), candidate_path.end(), [](Path &path1, Path &path2) {
            if(path1.size() == path2.size()) return path1 < path2;
            return path1.size() < path2.size();
        });

        Path chosen_path = candidate_path[0];
        reserve_path(chosen_path);
    }
}

void Greedy::BFS(int src, int dst) {
    vector<bool> visited(graph.get_num_nodes(), false);
    par.clear();
    par.resize(graph.get_num_nodes());
    for(int i = 0; i < graph.get_num_nodes(); i++) {
        par[i] = -1;
    }

    queue<int> que;
    que.push(src);
    visited[src] = true;

    if(graph.get_node_memory(src) < 1) return;
    while(!que.empty()) {
        int frt = que.front();
        que.pop();
        for(int neighbor : graph.adj_list[frt]) {
            int need_memory = (neighbor == dst) ? (1) : (2);
            if(graph.get_node_memory(neighbor) < need_memory) continue;
            if(!visited[neighbor]) {
                que.push(neighbor);
                visited[neighbor] = true;
                par[neighbor] = frt;
            }
        }
    }
}

Path Greedy::find_path(int src, int dst) {
    BFS(src, dst);
    if(par[dst] == -1) return {};
    Path path;
    int cur = dst;
    while(cur != -1) {
        path.push_back(cur);
        if(par[cur] != -1) assert(graph.adj_set[cur].count(par[cur]));
        cur = par[cur];
    }
    reverse(path.begin(), path.end());


    return path;
}