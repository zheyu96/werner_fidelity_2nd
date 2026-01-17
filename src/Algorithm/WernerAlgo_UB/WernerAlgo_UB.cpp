#include "WernerAlgo_UB.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

WernerAlgo_UB::WernerAlgo_UB(Graph graph, vector<pair<int,int>> requests, map<SDpair, vector<Path>> paths)
    : AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "WernerAlgo_UB";
}

void WernerAlgo_UB::variable_initialize() {
    int m = (int)requests.size() + graph.get_num_nodes() * graph.get_time_limit();

    double delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * m, 1.0 / epsilon));
    obj = m * delta;

    alpha.assign(requests.size(), delta);
    x.clear();
    x.resize(requests.size());

    int V = graph.get_num_nodes();
    int T_limit = graph.get_time_limit(); // 避免與 global T 混淆
    
    dpp.eps_bucket = graph.get_bucket_eps();
    double F_th = graph.get_fidelity_threshold();
    double w_th = (4.0 * F_th - 1.0) / 3.0;
    
    if(w_th <= 0) w_th = 1e-9;
    dpp.Zhat = sqrt(-log(w_th)) + 1e-9;
    dpp.Zmin = graph.get_Zmin();
    dpp.T    = time_limit - 1;
    dpp.tau_max = min(time_limit - 1, 5);
    dpp.eta  = graph.get_tao() / graph.get_time_limit();
    
    beta.assign(V, vector<double>(T_limit, INF));

    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T_limit; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? INF : (delta / cap);
        }
    }
}

Shape_vector WernerAlgo_UB::separation_oracle(){
    double most_violate = 1e9;
    Shape_vector todo_shape;
    
    for(int i = 0; i < (int)requests.size(); i++){
        int src = requests[i].first;
        int dst = requests[i].second;
        
        // [修正 1] 加上 const，解決 lvalue 錯誤
        const vector<Path>& req_paths = get_paths(src, dst); 
        
        for(const auto& path : req_paths){
            int T_search = dpp.T + 2; 
            int path_n = (int)path.size();
            
            DP_table.clear();
            DP_table.resize(T_search);
            for(int k = 0; k < T_search; k++){
                DP_table[k].resize(path_n);
                for(int j = 0; j < path_n; j++){
                    DP_table[k][j].resize(path_n);
                }
            }

            for(int t = 1; t <= (int)dpp.T; t++){
                run_dp_in_t(path, dpp, t);
                auto cur_val = eval_best_J(0, path_n - 1, t, alpha[i]);
                
                double path_prob = graph.path_Pr(path);
                if(cur_val.first != INF && (cur_val.first / path_prob < most_violate)){
                    most_violate = cur_val.first / path_prob;
                    todo_shape = backtrack_shape(cur_val.second, path);
                }
            }
        }
    }
    return todo_shape;
}

void WernerAlgo_UB::run_dp_in_t(const Path& path, const DPParam& dpp, int t) {
    const int n_len = (int)path.size();

    for(int a = 0; a < n_len - 1; a++) {
        for(int b = a + 1; b < n_len; b++) {
            int s = path[a], e = path[b];
            vector<ZLabel> cand;

            // 1. LEAF
            if(a + 1 == b) {
                for(int tlen = 1; tlen <= dpp.tau_max; tlen++) {
                    if (t - tlen < 0) continue;
                    double Bleaf = 0.0;
                    int st = t - tlen; 
                    for(int tt = st; tt <= t; tt++)
                        Bleaf += beta[s][tt] + beta[e][tt];
                    
                    double Zcur = 1.0L - (1.0L - graph.get_link_werner(s, e)) / tlen;
                    if (Zcur <= 0) Zcur = 1e-9;
                    double Zleaf = sqrt(-log(Zcur));
                    
                    if(Zleaf <= dpp.Zhat) {
                        ZLabel L(Bleaf, Zleaf, Op::LEAF, a, b, t, -1);
                        L.ent_time = {st, t};
                        cand.push_back(L);
                    }
                }
            }

            // 2. CONTINUE
            if(t > 0) {
                const auto& pre = DP_table[t-1][a][b];
                for(size_t p_id = 0; p_id < pre.size(); p_id++) {
                    double Zp = pre[p_id].Z + dpp.eta;
                    if(Zp <= dpp.Zhat) {
                        double Bp = pre[p_id].B + beta[s][t] + beta[e][t];
                        ZLabel L(Bp, Zp, Op::CONT, a, b, t, -1, (int)p_id);
                        cand.push_back(L);
                    }
                }
            }

            // 3. MERGE
            if(t > 0) {
                for(int k = a + 1; k < b; k++) {
                    const auto& L1 = DP_table[t-1][a][k];
                    const auto& L2 = DP_table[t-1][k][b];
                    
                    if(L1.empty() || L2.empty()) continue;
                    
                    for(size_t lid = 0; lid < L1.size(); lid++) {
                        for(size_t rid = 0; rid < L2.size(); rid++) {
                            const auto& left_seg = L1[lid];
                            const auto& right_seg = L2[rid];
                            
                            double Zp = sqrt((left_seg.Z + dpp.eta) * (left_seg.Z + dpp.eta) +
                                             (right_seg.Z + dpp.eta) * (right_seg.Z + dpp.eta));
                            
                            if(Zp <= dpp.Zhat) {
                                double Bp = left_seg.B + right_seg.B + beta[s][t] + beta[e][t];
                                ZLabel L(Bp, Zp, Op::MERGE, a, b, t, k, -1, (int)lid, (int)rid);
                                cand.push_back(L);
                            }
                        }
                    }
                }
            }
            bucket_by_Z(cand);
            DP_table[t][a][b] = cand;
        }
    }
}

void WernerAlgo_UB::pareto_prune_byZ(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    sort(cand.begin(), cand.end(), [](const ZLabel& x, const ZLabel& y){
        if(abs(x.Z - y.Z) > 1e-9) return x.Z < y.Z;
        return x.B < y.B;
    });

    vector<ZLabel> kept;
    double bestB = INF;
    for (auto& L : cand) {
        if (L.B + 1e-12 < bestB) {
            kept.push_back(L);
            bestB = L.B;
        }
    }
    cand.swap(kept);
}

void WernerAlgo_UB::bucket_by_Z(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    
    double q = 1 + dpp.eps_bucket;
    double invLogQ = 1.0 / log(q);
    map<int, ZLabel> buckets;

    for(auto& L : cand) {
        int k_idx;
        if(L.Z <= dpp.Zmin) k_idx = 0;
        else {
            double val = log(L.Z / dpp.Zmin) * invLogQ;
            k_idx = (int)floor(val + 1e-12);
            if(k_idx < 0) k_idx = 0;
        } 
        if(buckets.find(k_idx) == buckets.end() || L.B + 1e-12 < buckets[k_idx].B)
            buckets[k_idx] = L;
    }
    cand.clear();
    for(auto& P : buckets)
        cand.push_back(P.second);
    pareto_prune_byZ(cand);
}

Shape_vector WernerAlgo_UB::backtrack_shape(ZLabel leaf, const vector<int>& path){
    int left_id = path[leaf.a], right_id = path[leaf.b];
    
    if(leaf.op == Op::LEAF){
        Shape_vector result;
        if (leaf.ent_time.size() < 2) return Shape_vector{}; 
        result.push_back({left_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        result.push_back({right_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        return result;
    }
    
    if(leaf.op == Op::CONT){
        if (leaf.parent_id < 0 || leaf.parent_id >= (int)DP_table[leaf.t-1][leaf.a][leaf.b].size()) {
            return Shape_vector{};
        }
        ZLabel pre_label = DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
        Shape_vector last_time = backtrack_shape(pre_label, path);
        if (last_time.empty()) return Shape_vector{};
        last_time.front().second[0].second++;
        last_time.back().second[0].second++;
        return last_time;
    }
    
    if(leaf.op == Op::MERGE){
        if (leaf.k < 0) return Shape_vector{}; 
        ZLabel left_leaf = DP_table[leaf.t-1][leaf.a][leaf.k][leaf.left_id];
        Shape_vector left_result = backtrack_shape(left_leaf, path);
        ZLabel right_leaf = DP_table[leaf.t-1][leaf.k][leaf.b][leaf.right_id];
        Shape_vector right_result = backtrack_shape(right_leaf, path);

        if(left_result.empty() || right_result.empty()) return Shape_vector{};

        Shape_vector result = left_result;
        result.back().second.push_back(right_result.front().second.front());
        for(size_t i = 1; i < right_result.size(); i++) {
            result.push_back(right_result[i]);
        }
        result.front().second[0].second++;
        result.back().second[0].second++;
        return result;
    }
    return Shape_vector{};
}

int WernerAlgo_UB::split_dis(int s, int d, const ZLabel& L){
    if(L.op != Op::MERGE || L.k < 0) return 1000000000;
    int mid = (s + d) / 2;
    return abs(mid - L.k);
}

pair<double, WernerAlgo_UB::ZLabel> WernerAlgo_UB::eval_best_J(int s, int d, int t, double alp){
    double bestJ = 1e18;
    int bestdis = 1000000000;
    bool found = false;
    ZLabel tmp = {};

    for(const auto& L : DP_table[t][s][d]){
        double J = (alp + L.B) * exp(L.Z * L.Z);
        int dis = split_dis(s, d, L);
        
        if(J + EPS < bestJ || (fabs(J - bestJ) <= EPS && dis < bestdis)){
            bestJ = J;
            tmp = L;
            bestdis = dis;
            found = true;
        }
    }
    if(found) return {bestJ, tmp};
    else return {INF, tmp};
}

void WernerAlgo_UB::run() {
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        double eps = 1e-4;
        
        while (obj + eps < 1.0) {
            Shape_vector shape = separation_oracle();
            if (shape.empty()) break;

            double q = 1.0;
            for(int i = 0; i < (int)shape.size(); i++){
                map<int, int> need_amount;
                for(auto usedtime : shape[i].second){
                    for(int t = usedtime.first; t <= usedtime.second; t++)
                        need_amount[t]++;
                }
                for(auto P : need_amount){
                    int t = P.first;
                    double theta = P.second;
                    double cap = graph.get_node_memory_at(shape[i].first, t);
                    if (theta > 0)
                        q = min(q, cap / theta);
                }
            }

            if(q <= 1e-10) break;

            int req_idx = -1;
            int ln = shape.front().first;
            int rn = shape.back().first;
            
            for(int i = 0; i < (int)requests.size(); i++){
                if(requests[i] == make_pair(ln, rn)){
                    if(req_idx == -1 || alpha[i] > alpha[req_idx]){
                        req_idx = i;
                    }
                }
            }
            if(req_idx == -1) break;

            x[req_idx][shape] += q;

            double ori_alpha = alpha[req_idx];
            alpha[req_idx] = alpha[req_idx] * (1 + epsilon * q);
            obj += (alpha[req_idx] - ori_alpha);

            for(int i = 0; i < (int)shape.size(); i++){
                map<int, int> need_amount;
                for(auto usedtime : shape[i].second){
                    for(int t = usedtime.first; t <= usedtime.second; t++)
                        need_amount[t]++;
                }

                for(auto P : need_amount) {
                    int t = P.first;
                    int node_id = shape[i].first;
                    double theta = P.second;
                    double original_beta = beta[node_id][t];
                    double cap = graph.get_node_memory_at(node_id, t);
                    
                    if(cap > 0) {
                        beta[node_id][t] = beta[node_id][t] * (1 + epsilon * (q / (cap / theta)));
                        obj += (beta[node_id][t] - original_beta) * theta;
                    } else {
                        beta[node_id][t] = INF;
                    }
                }
            }
        } 

        double max_xim_sum = 0;
        double usage = 0;
        int memory_total_LP = 0;
        vector<bool> passed_node(graph.get_num_nodes(), false);

        res.clear();
        res["succ_request_cnt"] = 0;
        res["fidelity_gain"] = 0;
        res["pure_fidelity"] = 0;

        for(int i = 0; i < (int)requests.size(); i++) {
            double xim_sum = 0;
            for(auto P : x[i]) {
                double flow = P.second;
                xim_sum += flow;
                
                Shape shape_obj(P.first);
                double pr = graph.path_Pr(shape_obj);
                
                // [修正 2 & 3] 使用正確參數呼叫 get_fidelity
                // A, B, n, T, tao 應為全域變數 (config.h)
                double fid = shape_obj.get_fidelity(A, B, n, T, tao, graph.get_F_init());
                fid = ((1.0 + fid * 9.0) / 10.0);
                
                if(fid + EPS > graph.get_fidelity_threshold()) {
                    res["succ_request_cnt"] += flow * pr;
                    res["fidelity_gain"] += flow * (fid * pr);
                    res["pure_fidelity"] += flow * (fid * pr); 
                }

                for(auto id_mem : P.first) {
                    int node = id_mem.first;
                    if(!passed_node[node]) {
                        memory_total_LP += graph.get_node_memory(node);
                        passed_node[node] = true;
                    }
                    for(pair<int, int> mem_range : id_mem.second) {
                        usage += (mem_range.second - mem_range.first) * flow;
                    }
                }
            }
            max_xim_sum = max(max_xim_sum, xim_sum);
        }

        if(max_xim_sum > 0) {
            double total_capacity = (double)memory_total_LP * (double)graph.get_time_limit();
            if (total_capacity > 0)
                res["utilization"] = (usage / total_capacity) / max_xim_sum;
        }
        
        cerr << "\n[WernerAlgo_UB] UB Calculation Finished." << endl;
        cerr << "  Succ Req (UB): " << res["succ_request_cnt"] << endl;
        cerr << "  Fidelity (UB): " << res["fidelity_gain"] << endl;
    }
    cerr << "[" << algorithm_name << "] end" << endl;
}