#include "QCAST.h"

QCAST::QCAST() {
    method_name = "QCAST";
}

QCAST::~QCAST() {}

double QCAST::EXT(Path path_nodes, int w){
    vector<double> path_probs;
    for(int i = 1; i < (int)path_nodes.size(); i++) {
        path_probs.push_back(graph.get_entangle_succ_prob(path_nodes[i - 1], path_nodes[i]));
    }
    vector<vector<double>> Q(path_probs.size()); //Q[i][j] the prob of exactly j links success in ith edge on path 
    vector<vector<double>> P(path_probs.size());
    for(int i = 0; i < (int)path_probs.size(); i++){
        Q[i].resize(w + 1);
        P[i].resize(w + 1);
    }
    for(int k = 0;k < (int)path_probs.size();k++){
        for(int i=0;i<=w;i++){
            Q[k][i] = C(w, i) * pow(path_probs[k], i) * pow(1 - path_probs[k], w-i);
        }
    }
    for(int i=0;i<=w;i++){      //P0i = q0i
        P[0][i] = Q[0][i];
    }
    for(int k = 1; k < (int)path_probs.size(); k++){
        for(int i = 0; i <= w; i++){
            P[k][i] = 0;
            for(int l = i; l <= w; l++){
                P[k][i] += P[k-1][i] * Q[k][l];
            }
            for(int l = i + 1; l <= w; l++){
                P[k][i] += Q[k][i] * P[k-1][l];
            }
        }
    }

    //Et = q^h * sum(i * P[h][i])
    double ext = 0;
    for(int i = 1; i <= w; i++){
        ext += i * P.back()[i];
    }

    for(auto node : path_nodes) {
        ext *= graph.get_node_swap_prob(node);
    }
    return ext;
}

double QCAST::C(int n, int m){
    if(n < m){
        cerr<<"err:\tint QCAST combination m > n"<<endl;
        exit(1);
    }
    if(n == m || m == 0){return 1;}
    if(combination.find(make_pair(n, m)) != combination.end()){
        return combination[make_pair(n, m)];
    }
    combination[make_pair(n, m)] = C(n-1, m) + C(n-1, m-1);
    return combination[make_pair(n, m)];
}

const int maximum_major_path_per_request = 200;
const int maximum_path_length = 200;
const int maximum_total_number_of_path = 200;

pair<double, Path> QCAST::find_best_path(int src, int dst) {
    if(graph.get_node_memory(src) < 1 || paths[{src, dst}].size() > maximum_major_path_per_request){
        return {0, {}};
    }
    vector<double> dis(graph.get_num_nodes(), -1);
    vector<int> par(graph.get_num_nodes(), -1);

    priority_queue<pair<double, int>> pq;
    pq.emplace(INF, src);
    dis[src] = INF;
    while(!pq.empty()){
        int cur_node = pq.top().second;
        double distance = pq.top().first;
        pq.pop();
        if(dis[cur_node] > distance) continue;
        
        vector<int> path_nodes;
        int now = cur_node;
        while(now != -1){
            // cout<<"in find path: now = "<<now<<endl;
            path_nodes.push_back(now); 
            now = par[now];
        }
        reverse(path_nodes.begin(), path_nodes.end());
        if(path_nodes.size() > maximum_path_length){
            continue;
        }

        for(int neighbor : graph.adj_list[cur_node]){
            //cout<<ele<<" ";
            int need_memory = (neighbor == dst) ? (1) : (2);
            if(graph.get_node_memory(neighbor) < need_memory) continue;
            path_nodes.emplace_back(neighbor);
            int width = INF;
            for(int node : path_nodes) {
                width = min(width, graph.get_node_memory(node));
            }
            double ext = EXT(path_nodes, width);
            if(dis[neighbor] < ext){
                dis[neighbor] = ext;
                par[neighbor] = cur_node;
                pq.emplace(ext, neighbor);
            }
            path_nodes.pop_back();
        }
    }

    if(par[dst] == -1){
        //path is not found
        return {0, {}};
    }

    vector<int> path_nodes;
    vector<double> path_prob;
    int now = dst;
    while(now != -1){
        path_nodes.push_back(now); 

        now = par[now];
    }
    reverse(path_nodes.begin(), path_nodes.end());
    return {dis[dst], path_nodes};
}
void QCAST::build_paths(Graph _graph, vector<SDpair> _requests){
    if(DEBUG) cerr<< "---------QCAST::path_assignment----------" << endl;
    paths.clear();
    graph = _graph;
    requests = _requests;


    int total = 0;
    while(total < maximum_total_number_of_path){

        double mx_EXT = -1;
        Path mx_path;
        for(auto [src, dst] : requests){   //find the best path for every request
            double EXT;
            Path path;
            tie(EXT, path) = find_best_path(src, dst);
            if(EXT > mx_EXT) {
                mx_EXT = EXT;
                mx_path = path;
            }
        }
        //find the best path in requests
        if(mx_path.empty()) break;
        total++;
        reserve_path(mx_path);
    }

    if(DEBUG) cerr<< "---------QCAST::path_assignment----------end" << endl;
}