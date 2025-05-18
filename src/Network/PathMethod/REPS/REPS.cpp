#include "REPS.h"

using namespace std;

REPS::REPS() {
    method_name = "REPS";
}

REPS::~REPS() {}

void REPS::PFT_LP(vector<double> &t_plum, vector<map<pair<int, int>, double>> &f_plum) {

    //return value is store in t_plum and f_plum
    t_plum.clear();
    f_plum.clear();
    
    //do LP
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("OutputFlag", "0");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables
        
        
        vector<map<pair<int, int>, GRBVar>> f(requests.size());    //fi(u, v)
        for(int i = 0;i < (int)requests.size(); i++){
            for(int u = 0; u < graph.get_num_nodes(); u++){
                for(int v : graph.adj_list[u]){
                    f[i][make_pair(u, v)] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "f" + to_string(i) + "(" + to_string(u) + ", " + to_string(v) + ")");
                }
            }
        }


        vector<GRBVar> t;   //ti
        for(int i = 0; i < (int)requests.size(); i++){
            t.emplace_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t" + to_string(i)));
        }

        map<pair<int, int>, GRBVar> x; //x(u, v)
        for(int u = 0; u < graph.get_num_nodes(); u++){
            for(int v : graph.adj_list[u]) {
                x[make_pair(u, v)] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x(" + to_string(u) + ", " + to_string(v) + ")");
            }
        }

        // Set objective: 
        // maximize sum(ti)
        GRBLinExpr expr = 0;
        for (auto ti : t)
            expr += ti;
        
        model.setObjective(expr, GRB_MAXIMIZE);

        // Add constraint: 1(a) 
        for(int i=0;i<(int)requests.size();i++){
            GRBLinExpr expr = 0;
            int u = requests[i].first;
            for(int v : graph.adj_list[u]){
                expr += f[i][make_pair(u, v)];
                expr -= f[i][make_pair(v, u)];
            }
            model.addConstr(expr == t[i], "1a_" + to_string(i));
        }

        // Add constraint: 1(b) 
        for(int i=0;i<(int)requests.size();i++){
            expr = 0;
            int u = requests[i].second;
            for(int v : graph.adj_list[u]){
                expr += f[i][make_pair(u, v)];
                expr -= f[i][make_pair(v, u)];
            }
            model.addConstr(expr == GRBLinExpr(t[i], -1.0), "1b_"+ to_string(i));
        }

        // Add constraint: 1(c) 
        for(int i = 0; i < (int)requests.size(); i++){
            expr = 0;
            for(int u = 0; u < graph.get_num_nodes();u++){
                if(u == requests[i].first)continue;
                if(u == requests[i].second)continue;
                for(int v : graph.adj_list[u]){
                    expr += f[i][make_pair(u, v)];
                    expr -= f[i][make_pair(v, u)];
                }
                model.addConstr(expr == 0, "1c(" + to_string(i) + ", " + to_string(u) + ")");
            }
        }

        // Add constraint: 1(d) 
        for(int u = 0; u < graph.get_num_nodes(); u++){
            for(int v : graph.adj_list[u]){
                expr = 0;
                for(int i = 0; i < (int)requests.size(); i++){
                    expr += f[i][make_pair(u, v)];
                    expr += f[i][make_pair(v, u)];
                }
                double p = graph.get_entangle_succ_prob(u, v);
                model.addConstr(expr <= (x[make_pair(u, v)] * p), "1d(" + to_string(u) + ", " + to_string(v) + ")");
            }
        }

        // Add constraint: 1(e) 
        // for(int u = 0; u < graph.get_num_nodes(); u++){
        //     for(int v : graph.adj_list[u]){
        //         int c = INF;
        //         model.addConstr(x[make_pair(u, v)] <= c, "1e(" + to_string(u) + ", " + to_string(v) + ")");
        //     }
        // }

        // Add constraint: 1(f) 
        for(int  u = 0; u < graph.get_num_nodes(); u++){
            expr = 0;
            int m = graph.get_node_memory(u);
            for(int v : graph.adj_list[u]){
                expr += x[make_pair(u, v)];
            }
            model.addConstr(expr <= m, "1f(" + to_string(u) + ")");
        }

        // Optimize model
        model.optimize();

        //get t
        for(int i=0;i<(int)requests.size();i++){
            t_plum.emplace_back(t[i].get(GRB_DoubleAttr_X));
        }

        //get fi
        f_plum.resize(requests.size());
        for(int i = 0; i < (int)requests.size(); i++) {
            for(int u = 0; u < graph.get_num_nodes(); u++) {
                for(int v : graph.adj_list[u]) {
                    f_plum[i][make_pair(u, v)] = f[i][make_pair(u, v)].get(GRB_DoubleAttr_X);
                }
            }
        }

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
}

void REPS::build_paths(Graph _graph, vector<SDpair> _requests) {
    paths.clear();
    graph = _graph;
    requests = _requests;

    //PFT Using Progressive Rounding
    vector<double> t_plum;
    vector<map<pair<int, int>, double>> f_plum;
    vector<int> neighbor;
    bool flag = true;
    double width;
    vector<int> path;

    while(flag){
        PFT_LP(t_plum, f_plum);
        paths.clear();
        flag = false;

        vector<pair<int, Path>> tempPaths;
        for(int i = 0; i < (int)requests.size(); i++){
            int src = requests[i].first, dst = requests[i].second;
            while(true) {
                tie(path, width) = dijkstra(src, dst, f_plum[i]);
                if(path.empty() || width <= EPS) break;
                update_f_plum(path, width, f_plum[i]);
                tempPaths.emplace_back(width, path);
            }
        }

        for(auto [amount, path] : tempPaths) {
            reserve_path(path, amount);
        }
        // cout << "call PFT_LP in REPS::path_assignment()" << endl;
        // cout << "call PFT_LP in REPS::path_assignment()--end" << endl;
    }
}

pair<Path, double> REPS::dijkstra(int src, int dst, map<pair<int, int>, double> &f_plum_i) {
    vector<bool> vis(graph.get_num_nodes(), false);
    vector<double> flow(graph.get_num_nodes(), 0);
    vector<int> par(graph.get_num_nodes(), -1);
    priority_queue<pair<double, int>> pq;
    pq.push({INF, src});
    flow[src] = INF;
    while(!pq.empty()) {
        int cur = pq.top().second;
        pq.pop();
        if(vis[cur]) continue;
        vis[cur] = true;
        for(int v : graph.adj_list[cur]) {
            if(f_plum_i[{cur, v}] <= EPS) continue;
            double new_flow = min(flow[cur], f_plum_i[{cur, v}]);
            if(new_flow <= EPS) continue;
            if(flow[v] < new_flow) {
                flow[v] = new_flow;
                pq.push({flow[v], v});
                par[v] = cur;
            }
        }
    }

    if(!vis[dst]) return {{}, -1};

    Path path;
    int cur = dst;

    while(cur != -1) {
        path.push_back(cur);
        cur = par[cur];
    }

    reverse(path.begin(), path.end());
    return {path, flow[dst]};
}

void REPS::update_f_plum(Path path, double width, map<pair<int, int>, double>&f_plum_i) {
    for(int i = 1; i < (int)path.size(); i++) {
        int u = path[i - 1], v = path[i];
        if(f_plum_i[{u, v}] + EPS < width) {
            cerr << "error: f_plum_i is not enough" << endl;
            cerr << "f_plum = " << f_plum_i[{u, v}] << " width = " << width << endl;
            exit(1);
        }
        f_plum_i[{u, v}] -= width;
    }
}