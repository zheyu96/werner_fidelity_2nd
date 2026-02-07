#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <chrono> 
#include <omp.h>   
#include <cassert>

#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/WernerAlgo2_time/WernerAlgo2_time.h" 
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
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            int dist = graph.distance(i, j);
            if(dist >= length_lower && dist <= length_upper) {
                cand.emplace_back(i, j);
            }
        }
    }
    // 使用預設隨機引擎打亂
    shuffle(cand.begin(), cand.end(), generator);
    
    vector<SDpair> requests;
    uniform_int_distribution<int> unif(0, 1e9);
    for(SDpair sdpair : cand) {
        int cnt = unif(generator) % 4 + 3;
        while(cnt-- && (int)requests.size() < requests_cnt) requests.push_back(sdpair);
        if((int)requests.size() >= requests_cnt) break;
    }
    while((int)requests.size() < requests_cnt) requests.emplace_back(generate_new_request(n));
    return requests;
}

int main() {
    string file_path = "../data/";
    system("mkdir -p ../data/ans");

    // 1. 設定預設參數
    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["time_limit"] = 13;
    default_setting["avg_memory"] = 10;
    default_setting["tao"] = 0.002;
    default_setting["path_length"] = 4;
    default_setting["fidelity_threshold"] = 0.7;
    default_setting["epsilon"] = 0.75; 
    default_setting["Zmin"] = 0.02702867239;
    default_setting["bucket_eps"] = 0.01;
    default_setting["time_eta"] = 0.001;
    
    // 2. 設定要變動的參數範圍
    map<string, vector<double>> change_parameter;
    change_parameter["epsilon"] = {0.35, 0.45, 0.55, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}; 
    change_parameter["bucket_eps"] = {0.01, 0.02, 0.05, 0.1, 0.2};

    vector<string> X_names = {"epsilon", "bucket_eps"}; 
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt", "runtime"};
    int round = 1; 

    // 3. 預先產生固定請求
    vector<vector<SDpair>> default_requests(round);
    cerr << "Pre-generating requests..." << endl;
    for(int r = 0; r < round; r++) {
        string filename = file_path + "input/round_" + to_string(r) + ".input";
        Graph temp_g(filename, 20, 0.9, 10, 0.5, 0.98, 0.7, 0.25, 0.75, 2, 10, 0.002, 0.027, 0.01, 0.001);
        default_requests[r] = generate_requests(temp_g, (int)default_setting["request_cnt"], 3, 5);
    }

    // 4. 開始實驗 (多變數循環)
    for (const string& X_name : X_names) {
        if (change_parameter.find(X_name) == change_parameter.end()) continue;

        for(double change_value : change_parameter[X_name]) {
            cerr << "Testing " << X_name << " = " << change_value << endl;
            
            // 用來儲存每一輪的結果，外層 vector 對應 round
            // 內層 map 對應演算法名稱，再內層對應指標
            vector<map<string, map<string, double>>> result(round);

            #pragma omp parallel for
            for(int r = 0; r < round; r++) {
                string filename = file_path + "input/round_" + to_string(r) + ".input";
                
                // 複製預設參數並覆蓋當前變數
                map<string, double> input_param = default_setting;
                input_param[X_name] = change_value;

                Graph graph(filename, 
                            (int)input_param["time_limit"], 0.9, 
                            (int)input_param["avg_memory"], 0.5, 0.98, 
                            input_param["fidelity_threshold"], 0.25, 0.75, 2, 10, 
                            input_param["tao"], input_param["Zmin"], 
                            input_param["bucket_eps"], input_param["time_eta"]);

                vector<pair<int, int>> requests = default_requests[r];

                Greedy path_method;
                path_method.build_paths(graph, requests);
                auto paths = path_method.get_paths();

                // 建立演算法，確保傳入的是變動後的 epsilon
                AlgorithmBase* algo = new WernerAlgo2_time(graph, requests, paths, input_param["epsilon"],input_param["bucket_eps"]);

                auto start_time = chrono::high_resolution_clock::now();
                algo->run();
                auto end_time = chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed = end_time - start_time;

                string algo_name = algo->get_name(); 
                #pragma omp critical
                {
                    result[r][algo_name]["fidelity_gain"] = algo->get_res("fidelity_gain");
                    result[r][algo_name]["succ_request_cnt"] = algo->get_res("succ_request_cnt");
                    result[r][algo_name]["runtime"] = elapsed.count();
                }
                delete algo;
            }

            // 5. 寫入檔案
            for(string Y_name : Y_names) {
                string ans_filename = "ans/Greedy_" + X_name + "_" + Y_name + "_ZFA2.ans";
                ofstream ofs(file_path + ans_filename, ios::app);
                
                double sum = 0;
                string target_algo = "ZFA2_time"; 
                for(int r = 0; r < round; r++) sum += result[r][target_algo][Y_name];
                
                ofs << change_value << " " << sum / round << endl; 
                ofs.close();
            }
        }
    }
    return 0;
}