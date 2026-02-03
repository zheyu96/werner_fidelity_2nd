#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <chrono> // 用於高精度計時
#include <omp.h>   
#include <cassert>

#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/WernerAlgo_time/WernerAlgo_time.h" 
#include "Network/PathMethod/Greedy/Greedy.h"

using namespace std;

// --- 輔助函式 ---
SDpair generate_new_request(int num_of_node){
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node-1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);
    return make_pair(node1, node2);
}

vector<SDpair> generate_requests(Graph graph, int requests_cnt, int length_lower, int length_upper) {
    int n = graph.get_num_nodes();
    vector<SDpair> cand;
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            int dist = graph.distance(i, j);
            if(dist >= length_lower && dist <= length_upper) {
                cand.emplace_back(i, j);
            }
        }
    }
    random_shuffle(cand.begin(), cand.end());
    vector<SDpair> requests;
    for(SDpair sdpair : cand) {
        int cnt = unif(generator) % 4 + 3;
        while(cnt--) requests.push_back(sdpair);
    }
    while((int)requests.size() < requests_cnt) requests.emplace_back(generate_new_request(n));
    while((int)requests.size() > requests_cnt) requests.pop_back();
    return requests;
}

// --- 主程式 ---
int main() {
    string file_path = "../data/";
    
    // 確保 ans 目錄存在（建議手動建立，或透過系統指令）
    system("mkdir -p ../data/ans");

    // 1. 設定預設參數
    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["time_limit"] = 20;
    default_setting["avg_memory"] = 10;
    default_setting["tao"] = 0.002;
    default_setting["path_length"] = 4;
    default_setting["fidelity_threshold"] = 0.7;
    default_setting["epsilon"] = 0.35; 
    default_setting["Zmin"] = 0.02702867239;
    default_setting["bucket_eps"] = 0.01;
    default_setting["time_eta"] = 0.001;
    
    // 2. 設定實驗變數：測試 epsilon 對 runtime 的影響
    map<string, vector<double>> change_parameter;
    change_parameter["epsilon"] = {0.35, 0.45, 0.55, 0.65, 0.75, 0.85,0.95}; 

    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt", "runtime"};
    int round = 5; 

    // 3. 預先產生固定請求，確保實驗公平性
    vector<vector<SDpair>> default_requests(round);
    cerr << "Pre-generating requests..." << endl;
    for(int r = 0; r < round; r++) {
        string filename = file_path + "input/round_" + to_string(r) + ".input";
        // 建立暫時的 Graph 來計算距離並生成請求
        Graph temp_g(filename, 20, 0.9, 10, 0.5, 0.98, 0.7, 0.25, 0.75, 2, 10, 0.002, 0.027, 0.01, 0.001);
        default_requests[r] = generate_requests(temp_g, (int)default_setting["request_cnt"], 3, 5);
    }

    // 4. 開始實驗
    string X_name = "epsilon"; 
    for(double change_value : change_parameter[X_name]) {
        cerr << "Testing " << X_name << " = " << change_value << endl;
        vector<map<string, map<string, double>>> result(round);

        #pragma omp parallel for
        for(int r = 0; r < round; r++) {
            string filename = file_path + "input/round_" + to_string(r) + ".input";
            map<string, double> input_parameter = default_setting;
            input_parameter[X_name] = change_value;

            // 初始化當前實驗的 Graph
            Graph graph(filename, 
                        (int)input_parameter["time_limit"], 0.9, 
                        (int)input_parameter["avg_memory"], 0.5, 0.98, 
                        input_parameter["fidelity_threshold"], 0.25, 0.75, 2, 10, 
                        input_parameter["tao"], input_parameter["Zmin"], 
                        input_parameter["bucket_eps"], input_parameter["time_eta"]);

            // 使用預先生成的 requests
            vector<pair<int, int>> requests = default_requests[r];

            // 尋找路徑
            Greedy path_method;
            path_method.build_paths(graph, requests);
            map<SDpair, vector<Path>> paths = path_method.get_paths();

            // 建立演算法並傳入測試的 epsilon
            double eps_to_use = change_value;
            AlgorithmBase* algo = new WernerAlgo_time(graph, requests, paths, eps_to_use);

            // 計時開始
            auto start_time = chrono::high_resolution_clock::now();
            algo->run();
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end_time - start_time;

            // 紀錄結果
            string name = algo->get_name(); // 預期為 "ZFA_time"
            #pragma omp critical
            {
                result[r][name]["fidelity_gain"] = algo->get_res("fidelity_gain");
                result[r][name]["succ_request_cnt"] = algo->get_res("succ_request_cnt");
                result[r][name]["runtime"] = elapsed.count();
            }
            delete algo;
        }

        // 5. 寫入平均結果
        for(string Y_name : Y_names) {
            string ans_filename = "ans/Greedy_" + X_name + "_" + Y_name + ".ans";
            ofstream ofs(file_path + ans_filename, ios::app);
            ofs << change_value << " ";
            
            string algo_name = "ZFA_time"; 
            double sum = 0;
            for(int r = 0; r < round; r++) sum += result[r][algo_name][Y_name];
            ofs << sum / round << endl; 
            ofs.close();
        }
    }
    return 0;
}