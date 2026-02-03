#include "WernerAlgo2_time.h"
#include <fstream> 
#include <iostream>
#include <cmath>

using namespace std;

WernerAlgo2_time::WernerAlgo2_time(Graph graph,vector<pair<int,int>> requests,map<SDpair, vector<Path>> paths,double eps): AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "ZFA2_time";
    this->epsilon=eps;
}

void WernerAlgo2_time::variable_initialize() {
    // 與 MyAlgo1 類似：初始化 dual 與目標
    int m = (int)requests.size()
          + graph.get_num_nodes() * graph.get_time_limit();

    double delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * m, 1.0 / epsilon));
    obj = m * delta;

    alpha.assign(requests.size(), delta);
    x.clear();
    x.resize(requests.size());
    int V = graph.get_num_nodes();
    int T = graph.get_time_limit();
    dpp.eps_bucket = graph.get_bucket_eps();
    double F_th=graph.get_fidelity_threshold();
    double w_th=(4.0*F_th-1.0)/3.0;
    dpp.Zhat = sqrt(-log(w_th))+1e-9;
    dpp.Zmin = graph.get_Zmin();
    dpp.T    = time_limit-1;
    dpp.tau_max=min(time_limit-1,5);
    dpp.eta  = graph.get_tao()/graph.get_time_limit();
    beta.assign(V, vector<double>(T, INF));

    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? INF : (delta / cap);
        }
    }
}

Shape_vector WernerAlgo2_time::separation_oracle(){
    // 針對每個 request 找最小成本 shape，選最好的回傳
    double most_violate=1e9;
    bool flag=0;
    Shape_vector todo_shape;
    for(int i=0;i<requests.size();i++){
        int src=requests[i].first,dst=requests[i].second;
        vector<Path> paths=get_paths(src,dst);
        for(int p=0;p<paths.size();p++){
            // Example initialization: adjust T and n as needed for your context
            int T=dpp.T+5;
            int n=paths[p].size()+5;
            DP_table.clear();
            DP_table.resize(T);
            for(int i=0;i<DP_table.size();i++){
                DP_table[i].resize(n);
                for(int j=0;j<DP_table[i].size();j++)
                    DP_table[i][j].resize(n);
            }
            for(int t=1;t<=dpp.T;t++){
                run_dp_in_t(paths[p],dpp,t);
                auto cur_val=eval_best_J(0,paths[p].size()-1,t,alpha[i]);
                if(cur_val.first/graph.path_Pr(paths[p])<most_violate){
                    most_violate=cur_val.first/graph.path_Pr(paths[p]);
                    todo_shape=backtrack_shape(cur_val.second,paths[p]);
                }
            }
        }
    }
    return todo_shape;
}

void WernerAlgo2_time::run_dp_in_t(const Path& path, const DPParam& dpp,int t) {
    const int T = graph.get_time_limit();
    const int n = (int)path.size();
    long long total_before = 0;
    long long total_after = 0;
    // -------- t = 1..T-1 外圈時間迴圈 --------
    for(int a=0;a<n-1;a++)
        for(int b=a+1;b<n;b++){
            int s=path[a],e=path[b];
            vector<ZLabel> cand;
            //leaf
            if(a+1==b){
                for(int tlen=1;tlen<=dpp.tau_max;tlen++){
                    if (t - tlen < 1) continue;
                    double Bleaf=0.0;
                    int st=t-tlen;
                    if(st<1)continue;
                    for(int tt=st;tt<=t;tt++)
                        Bleaf+=beta[s][tt]+beta[e][tt];
                    double Zcur=1.0L-(1.0L-graph.get_link_werner(s,e))/tlen;
                    double Zleaf=sqrt(-log(Zcur));
                    if(Zleaf<=dpp.Zhat){
                        ZLabel L(Bleaf,Zleaf,Op::LEAF,a,b,t,-1);
                        L.ent_time={t-tlen,t};
                        cand.push_back(L);
                    }
                }
            }
            //continue
            auto pre=DP_table[t-1][a][b];
            for(int p_id=0;p_id<pre.size();p_id++){
                double Zp=pre[p_id].Z+dpp.eta;
                if(Zp<=dpp.Zhat){
                    double Bp=pre[p_id].B+beta[s][t]+beta[e][t];
                    ZLabel L(Bp,Zp,Op::CONT,a,b,t,-1,p_id);
                    cand.push_back(L);
                }
            }
            //merge
            for(int k=a+1;k<b;k++){
                auto L1=DP_table[t-1][a][k],L2=DP_table[t-1][k][b];
                if(L1.size()==0||L2.size()==0) continue;
                for(int lid=0;lid<L1.size();lid++)
                    for(int rid=0;rid<L2.size();rid++){
                        auto left_seg=L1[lid],right_seg=L2[rid];
                        double Zp=sqrt((left_seg.Z+dpp.eta)*(left_seg.Z+dpp.eta)+
                                        (right_seg.Z+dpp.eta)*(right_seg.Z+dpp.eta));
                        if(Zp<=dpp.Zhat){
                            double Bp=left_seg.B+right_seg.B+beta[s][t]+beta[e][t];
                            ZLabel L(Bp,Zp,Op::MERGE,a,b,t,k,-1,lid,rid);
                            cand.push_back(L);
                        }
                    }
            }
            total_before = cand.size();
            bucket_by_Z(cand);
            total_after = cand.size();
            cout << "DP_table[" << t << "][" << a << "][" << b << "] - before: " << total_before << ", after: " << total_after << endl;
            cout << "Reduction ratio: " << (total_before == 0 ? 0.0 : (double)total_after / total_before) << endl;
            for(auto i:cand){
                cout << "Z: " << i.Z << ", B: " << i.B << endl;
            }
            // Convert cand (vector<ZLabel>) to vector<shared_ptr<ZLabel>>
            DP_table[t][a][b]=cand;
        }
}

void WernerAlgo2_time::pareto_prune_byZ(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    sort(cand.begin(), cand.end(), [](const ZLabel& x, const ZLabel& y){
        if(x.Z!=y.Z) return x.Z < y.Z;
        return x.B<y.B;
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

void WernerAlgo2_time::bucket_by_Z(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    double q=1+dpp.eps_bucket;
    double invLogQ=1.0/log(q);
    map<double,ZLabel> buckets;
    for(auto L:cand){
        double k;
        if(L.Z<=dpp.Zmin) k=0.0;
        else{
            k=floor(log(L.Z/dpp.Zmin)*invLogQ+1e-12);
            if(k<0) k=0.0;
        } 
        if(buckets.find(k)==buckets.end()||L.B+1e-12<buckets[k].B)
            buckets[k]=L;
    }
    vector<ZLabel> bucketed;
    for(auto L:buckets)
        bucketed.push_back(L.second);
    pareto_prune_byZ(bucketed);
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y){
        return x.Z < y.Z;
    });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo2_time::backtrack_shape(ZLabel leaf,const vector<int>& path){
    int left_id=path[leaf.a],right_id=path[leaf.b];
    if(leaf.op==Op::LEAF){
        Shape_vector result;
        if (leaf.ent_time.size() < 2) return Shape_vector{}; 
        assert(leaf.ent_time.size() == 2);
        // 這邊的 fidelity 計算僅為 backtrack 過程中的檢查，與主要統計邏輯分開
        double fid=graph.get_F_init(left_id,right_id);
        double w=graph.get_link_werner(left_id,right_id);
        double new_w=1.0L-(1.0L-graph.get_link_werner(left_id,right_id))/(leaf.ent_time[1]-leaf.ent_time[0]);
        double new_fid= (3.0*new_w +1.0)/4.0;
        
        result.push_back({left_id,{{leaf.ent_time[0],leaf.ent_time[1]}}});
        result.push_back({right_id,{{leaf.ent_time[0],leaf.ent_time[1]}}});
        return result;
    }
    if(leaf.op==Op::CONT){
        assert(leaf.parent_id>=0&&leaf.parent_id<DP_table[leaf.t-1][leaf.a][leaf.b].size());
        ZLabel pre_label=DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
        Shape_vector last_time=backtrack_shape(pre_label,path);
        if (last_time.empty()) return Shape_vector{};
        auto & prel=last_time.front().second[0],&prer=last_time.back().second[0];
        assert(last_time.front().first==path[leaf.a]);
        assert(last_time.back().first==path[leaf.b]);
        assert(prel.second==leaf.t-1);
        assert(prer.second==leaf.t-1);
        prel.second++;
        prer.second++;
        return last_time;
    }
    if(leaf.op==Op::MERGE){
        Shape_vector left_result,right_result,result;
        assert(leaf.k>=0);
        int k_id=path[leaf.k];
        ZLabel left_leaf=DP_table[leaf.t-1][leaf.a][leaf.k][leaf.left_id];
        left_result=backtrack_shape(left_leaf,path);
        ZLabel right_leaf=DP_table[leaf.t-1][leaf.k][leaf.b][leaf.right_id];
        right_result=backtrack_shape(right_leaf,path);
        if(DEBUG) {
            assert(left_result.front().first == path[leaf.a]);
            assert(left_result.front().second[0].second == leaf.t - 1);
            assert(left_result.front().second.size() == 1);
            assert(left_result.back().first == k_id);
            assert(right_result.front().first == k_id);
            assert(right_result.back().first == path[leaf.b]);
            assert(right_result.back().second[0].second == leaf.t - 1);
            assert(left_result.back().second.size() == 1);
        }

        for(int i = 0; i < (int)left_result.size(); i++) {
            result.push_back(left_result[i]);
        }
        result.back().second.push_back(right_result.front().second.front());
        for(int i = 1; i < (int)right_result.size(); i++) {
            result.push_back(right_result[i]);
        }

        result.front().second[0].second++;
        result.back().second[0].second++;
        return result;
    }
    return Shape_vector{};
    // Handle unexpected Op value
    cerr << "[WernerAlgo2_time::backtrack_shape] Warning: Unknown Op value encountered." << std::endl;
}
int WernerAlgo2_time::split_dis(int s,int d,WernerAlgo2_time::ZLabel& L){
    if(L.op!=WernerAlgo2_time::Op::MERGE||L.k<0) return 1e18/4;
    int mid=(s+d)/2;
    return abs(mid-L.k);
}
pair<double,WernerAlgo2_time::ZLabel> WernerAlgo2_time::eval_best_J(int s, int d, int t, double alp){
    double bestJ=1e18;
    int bestdis=1e18/4;
    int flag=0;
    ZLabel tmp={};
    for(auto L:DP_table[t][s][d]){
        double J=(alp+L.B)*exp(L.Z*L.Z);
        int dis=split_dis(s,d,L);
        if(J+EPS<bestJ||(fabs(J-bestJ)<=EPS&&dis<bestdis)){
            bestJ=J;
            tmp=L;
            bestdis=dis;
            flag=1;
        }
    }
    if(flag) return {bestJ,tmp};
    else return {INF,tmp};
}

void WernerAlgo2_time::run() {
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        //cerr << "\033[1;31m"<< "[WernerAlgo's parameter] : "<< dpp.Zmin<<" "<<dpp.eps_bucket<<" "<<dpp.eta<< "\033[0m"<< endl;
        int it=0;
        double eps=1e-4;
        while (obj+eps < 1.0) {
            //if(++it>200) break;
            Shape_vector shape=separation_oracle();
            if (shape.empty()) break;
            // 先用MyAlgo1的框架刻出來
            double q = 1.0;
            for(int i=0;i<shape.size();i++){
                map<int,int> need_amount;
                for(pair<int,int> usedtime:shape[i].second){
                    int start=usedtime.first,end=usedtime.second;
                    for(int t=start;t<=end;t++)
                        need_amount[t]++;
                }
                for(pair<int,int>P:need_amount){
                    int t=P.first;
                    double theta=P.second;
                    q=min(q,graph.get_node_memory_at(shape[i].first,t)/theta);
                }
            }
            if(q<=1e-10) break;
            int req_idx=-1;
            for(int i=0;i<requests.size();i++){
                int ln=shape.front().first,rn=shape.back().first;
                if(requests[i]==make_pair(ln,rn)){
                    if(req_idx==-1||alpha[req_idx]>alpha[i]){
                        req_idx=i;
                    }
                }
            }
            if(req_idx==-1) break;
            x[req_idx][shape]+=q;
            double ori=alpha[req_idx];
            alpha[req_idx]=alpha[req_idx]*(1+epsilon*q);
            obj+=(alpha[req_idx]-ori);
            for(int i=0;i<shape.size();i++){
                map<int,int> need_amount;
                for(pair<int,int> usedtime:shape[i].second){
                    int start=usedtime.first,end=usedtime.second;
                    for(int t=start;t<=end;t++)
                        need_amount[t]++;
                }

                for(pair<int, int> P : need_amount) {
                    int t = P.first;
                    int node_id = shape[i].first;
                    double theta = P.second;
                    double original = beta[node_id][t];
                    if(graph.get_node_memory_at(node_id, t) == 0) {
                        beta[node_id][t] = INF;
                    } else {
                        beta[node_id][t] = beta[node_id][t] * (1 + epsilon * (q / (graph.get_node_memory_at(node_id, t) / theta)));
                    }
                    obj += (beta[node_id][t] - original) * theta;
                }
            }
        }
        vector<pair<double, Shape_vector>> shapes;

        for(int i = 0; i < (int)requests.size(); i++) {
            for(auto P : x[i]) {
                shapes.push_back({P.second, P.first});
            }
        }

        sort(shapes.begin(), shapes.end(), [](pair<double, Shape_vector> left, pair<double, Shape_vector> right) {
            return left.first > right.first;
        });

        vector<bool> used(requests.size(), false);
        vector<int> finished;
        
        // [新增] 統計資料結構：Key=Duration, Value=List of {W_raw, W_new}
        map<int, vector<pair<double, double>>> purification_stats;

        for(pair<double, Shape_vector> P : shapes) {
            Shape shape = Shape(P.second);
            int request_index = -1;
            for(int i = 0; i < (int)requests.size(); i++) {
                if(used[i] == false && requests[i] == make_pair(shape.get_node_mem_range().front().first, shape.get_node_mem_range().back().first)) {
                    request_index = i;
                }
            }

            if(request_index == -1 || used[request_index]) continue;
            
            // 呼叫 check_resource 時開啟純化 (Enable Purification)
            if(graph.check_resource(shape, true, true)) {
                used[request_index] = true;
                graph.reserve_shape(shape, true);
                finished.push_back(request_index);

                // [新增] 收集純化資訊
                Shape_vector sv = shape.get_node_mem_range();
                // 遍歷每一段 Link (從 node i 到 node i+1)
                for(size_t i = 0; i < sv.size() - 1; ++i) {
                    int u = sv[i].first;
                    int v = sv[i+1].first;

                    // Leaf node 的 memory range 最後一個區間代表與下一個節點的 Entanglement
                    int start = sv[i].second.back().first;
                    int end = sv[i].second.back().second;
                    int duration = end - start;

                    // 只有當 duration > 1 代表有進行純化操作
                    if (duration > 1) {
                        double w_raw = graph.get_link_werner(u, v);
                        double w_new = 1.0 - (1.0 - w_raw) / (double)duration;
                        purification_stats[duration].push_back({w_raw, w_new});
                    }
                }
            }
        }

        sort(finished.rbegin(), finished.rend());
        for(auto fin : finished) {
            requests.erase(requests.begin() + fin);
        }

        // [新增] 將統計結果寫入檔案
        // 檔案路徑與 main.cpp 中的 log 路徑一致 (假設在 ../data/log/)
        // 使用 append 模式，避免覆蓋之前的回合記錄
        string log_file_path = "../data/log/ZFA2_Purification_Stats.txt";
        ofstream log_file(log_file_path, ios::app);
        
        if (log_file.is_open()) {
            if (!purification_stats.empty()) {
                log_file << "--- Round Log ---" << endl;
                for(auto const& [len, vec] : purification_stats) {
                    log_file << "Duration " << len << " (Count: " << vec.size() << "):" << endl;
                    for(auto const& pair : vec) {
                        double w_raw = pair.first;
                        double w_new = pair.second;
                        // 轉換為 Fidelity 方便閱讀: F = (3W+1)/4
                        double f_raw = (3.0 * w_raw + 1.0) / 4.0;
                        double f_new = (3.0 * w_new + 1.0) / 4.0;
                        
                        log_file << "  W: " << w_raw << " -> " << w_new 
                                 << " | F: " << f_raw << " -> " << f_new << endl;
                    }
                }
                log_file << "-----------------" << endl;
            }
            log_file.close();
        } else {
            cerr << "[Warning] Unable to open log file: " << log_file_path << endl;
        }

        // 同時印到 cerr 方便即時觀察
        cerr << "\n=== [WernerAlgo2 Purification Stats] ===" << endl;
        if (purification_stats.empty()) {
            cerr << "No links were purified (all duration = 1)." << endl;
        } else {
            for(auto const& [len, vec] : purification_stats) {
                cerr << "Duration " << len << " (Count: " << vec.size() << "):" << endl;
                // 為了避免洗版，終端機只印數量，詳細內容看檔案
            }
        }
        cerr << "========================================\n" << endl;
    }
    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}