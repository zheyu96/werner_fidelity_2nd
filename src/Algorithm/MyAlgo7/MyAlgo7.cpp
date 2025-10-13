#include "MyAlgo7.h"
MyAlgo7::MyAlgo7(Graph graph, vector<pair<int, int>> requests, map<SDpair, vector<Path>> paths):
    AlgorithmBase(graph, requests, paths) {
    algorithm_name = "MyAlgo7";

}
void variable_initialize(){
    epsilon=0.9;
    double m=request.size()+(double)graph.get_num_nodes()*(double)graph.get_time_limit();
}
Shape_vector MyAlgo1::separation_oracle() {
    Shape_vector min_shape;
    double min_value = INF;
    for(int i = 0; i < (int)requests.size(); i++) {
        int src = requests[i].first, dst = requests[i].second;
        // cerr << "[MyAlgo1] " << "path len = " << graph.get_path(src, dst).size() << endl;
        auto result = find_min_shape(src, dst, alpha[i]);
        Shape_vector shape = result.first;
        double value = result.second;
        if(value < min_value) {
            min_shape = shape;
            min_value = value;
        }
        // cerr << "[MyAlgo1] " << "find shape" << endl;
    }
    return min_shape;
}
strategy path_dp(double lanbda,int path_num,int time_len){
    strategy cur=new strategy(path_num,time_len);

}
pair<Shape_vector, double> MyAlgo7::find_min_shape(int src, int dst, double alp) {
    vector<Path> paths = get_paths(src, dst);
    Shape_vector best_shape;
    double best_cost = INF;
    for(Path path : paths) {
        double lambda_left=0,lambda_right=TBA;
        
        dp_table.clear(0);
        dp.clear();
        dp.resize(path.size());
        par.clear();
        par.resize(path.size());
        caled.clear();
        caled.resize(path.size());
        for(int i = 0; i < (int)path.size(); i++) {
            dp[i].resize(path.size());
            par[i].resize(path.size());
            caled[i].resize(path.size());
            for(int j = 0; j < (int)path.size(); j++) {
                dp[i][j].resize(time_limit,{INF,INF});
                par[i][j].resize(time_limit, -2);
                caled[i][j].resize(time_limit, false);
            }
        }

        double best = INF;
        int best_time = -1;
        for(int t = 0; t < time_limit; t++) {
            double result = recursion_calculate_min_shape(0, path.size() - 1, t, path,lambda);
            result = (result + alp) / graph.path_Pr(path);
            if(best > result) {
                best = result;
                best_time = t;
            }
        }

        if(best_time == -1) continue;

        if(best < best_cost) {
            best_shape = recursion_find_shape(0, (int)path.size() - 1, best_time, path);
            best_cost = best;
        }
    }

    if(best_cost == INF) return {{}, INF};
    return {best_shape, best_cost};
}
state MyAlgo7::recursion_calculate_min_shape(int left, int right, int t, vector<int> &path,double lambda) {
    if(t <= 0) return {INF,INF};
    // cerr << left << " " << right << " " << t << " " << (int)path.size() << endl;

    int left_id = path[left], right_id = path[right];
    //leaf
    if(left == right - 1) {
        state leaf;
        leaf.B=beta[left_id][t - 1] + beta[left_id][t] + beta[right_id][t - 1] + beta[right_id][t];
        leaf.W=graph.get_edge_W(left_id,right_id);
        return leaf;
    }

    //calculate before
    if(caled[left][right][t]) return dp[left][right][t];
    
    //continue
    state best = recursion_calculate_min_shape(left, right, t - 1, path);
    
    //merge
    int best_k = -1;
    for(int k = left + 1; k < right; k++) {
        state left_result=recursion_calculate_min_shape(left, k, t - 1, path);
        state right_result=recursion_calculate_min_shape(k, right, t - 1, path);
        state result=left_result + right_result;
        if(eval(result,lambda)+EPS<eval(best,lambda)){
            best=result;
            best_k=k;
        }
    }

    caled[left][right][t] = true;
    par[left][right][t] = best_k;
    return dp[left][right][t] = best;
}

void MyAlgo7::run(){
    int round=1;
    while(round--&&!request.empty()){
        variable_initialize();

    }
}