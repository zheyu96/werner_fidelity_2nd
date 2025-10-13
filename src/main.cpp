
#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/MyAlgo1/MyAlgo1.h"
#include "Algorithm/MyAlgo2/MyAlgo2.h"
#include "Algorithm/MyAlgo3/MyAlgo3.h"
#include "Algorithm/MyAlgo4/MyAlgo4.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
#include "Algorithm/MyAlgo6/MyAlgo6.h"
#include "Algorithm/WernerAlgo/WernerAlgo.h"
#include "Network/PathMethod/PathMethodBase/PathMethod.h"
#include "Network/PathMethod/Greedy/Greedy.h"
#include "Network/PathMethod/QCAST/QCAST.h"
#include "Network/PathMethod/REPS/REPS.h"

using namespace std;

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

    while((int)requests.size() < requests_cnt) {
        requests.emplace_back(generate_new_request(n));
    }

    while((int)requests.size() > requests_cnt) {
        requests.pop_back();
    }

    return requests;
}

int main(){
    string file_path = "../data/";

    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["entangle_lambda"] = 0.045;
    default_setting["time_limit"] = 13;
    default_setting["avg_memory"] = 10; // 16
    default_setting["tao"] = 0.002;
    default_setting["path_length"] = 5;
    default_setting["min_fidelity"] = 0.7;
    default_setting["max_fidelity"] = 0.98;
    default_setting["swap_prob"] = 0.9;
    default_setting["fidelity_threshold"] = 0.5;
    default_setting["entangle_time"] = 0.00025;
    default_setting["entangle_prob"] = 0.01;

    map<string, vector<double>> change_parameter;
    change_parameter["request_cnt"] = {10, 30, 50, 70, 90};
    change_parameter["num_nodes"] = {40, 70, 100, 130, 160};
    change_parameter["min_fidelity"] = {0.6, 0.7, 0.8, 0.9, 0.95};
    change_parameter["avg_memory"] = {6, 8, 10, 12, 14};
    // change_parameter["tao"] = {0.3, 0.4, 0.5, 0.6, 0.7};
    change_parameter["tao"] = {0.00125, 0.0015, 0.002, 0.00225, 0.0023};
    change_parameter["path_length"] = {3, 6, 9, 12, 15};
    change_parameter["swap_prob"] = {0.5, 0.6, 0.7, 0.8, 0.9};
    change_parameter["fidelity_threshold"] = {0.4, 0.45, 0.5, 0.55, 0.6};
    change_parameter["time_limit"] = {5, 9, 13, 17, 21};
    change_parameter["entangle_lambda"] = {0.025, 0.035, 0.045, 0.055, 0.065};
    change_parameter["entangle_time"] = {0.0001, 0.00025, 0.0004, 0.00055, 0.0007};
    change_parameter["entangle_prob"] = {0.0001, 0.001, 0.01, 0.1, 1};

    int round = 50;
    vector<vector<SDpair>> default_requests(round);

    #pragma omp parallel for
    for(int r = 0; r < round; r++) {
        int num_nodes = default_setting["num_nodes"];
        int avg_memory = default_setting["avg_memory"];
        // int request_cnt = default_setting["request_cnt"];
        int time_limit = default_setting["time_limit"];
        double min_fidelity = default_setting["min_fidelity"];
        double max_fidelity = default_setting["max_fidelity"];

        double swap_prob = default_setting["swap_prob"];
        double fidelity_threshold = default_setting["fidelity_threshold"];
        int length_upper = default_setting["path_length"] + 1;
        int length_lower = default_setting["path_length"] - 1;
        for(string X_name : {"entangle_prob"}) {
            map<string, double> input_parameter = default_setting;
            for(double change_value : change_parameter[X_name]) {
                vector<map<string, map<string, double>>> result(round);
                input_parameter[X_name] = change_value;
                // double entangle_lambda = input_parameter["entangle_lambda"];
                // double entangle_time = input_parameter["entangle_time"];
                double entangle_prob = input_parameter["entangle_prob"];
                string filename = file_path + "input/round_" + to_string(r) + ".input";
                string command = "python3 graph_generator.py ";
                double A = 0.25, B = 0.75, tao = default_setting["tao"], T = 10, n = 2;
                // derandom
                string parameter = to_string(num_nodes);
                cerr << (command + filename + " " + parameter) << endl;
                if(system((command + filename + " " + parameter).c_str()) != 0){
                    cerr<<"error:\tsystem proccess python error"<<endl;
                    exit(1);
                }
                Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao, entangle_prob);
                default_requests[r] = generate_requests(graph, 100, length_lower, length_upper);
            }
        }
    }




    // vector<string> X_names = {"time_limit", "request_cnt", "num_nodes", "avg_memory", "tao"};
    vector<string> X_names = {"entangle_time", "request_cnt", "time_limit", "tao", "fidelity_threshold", "avg_memory", "min_fidelity", "entangle_lambda", "swap_prob"};
    //vector<string> X_names = {"entangle_prob","request_cnt"};
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt"};
    vector<string> algo_names = {"WernerAlgo","MyAlgo1", "MyAlgo2", "MyAlgo3", "Merge", "Linear", "ASAP"};
    // init result


    vector<PathMethod*> path_methods;
    path_methods.emplace_back(new Greedy());
    path_methods.emplace_back(new QCAST());
    path_methods.emplace_back(new REPS());
    for(PathMethod *path_method : path_methods) {

        for(string X_name : X_names) {
            for(string Y_name : Y_names){
                string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                fstream file( file_path + filename, ios::out );
            }
        }

        for(string X_name : X_names) {
            /* if(path_method->get_name() != "Greedy" && X_name != "request_cnt")
                continue; */
                
            map<string, double> input_parameter = default_setting;

            for(double change_value : change_parameter[X_name]) {
                vector<map<string, map<string, double>>> result(round);
                input_parameter[X_name] = change_value;

                // int num_nodes = input_parameter["num_nodes"];
                int avg_memory = input_parameter["avg_memory"];
                int request_cnt = input_parameter["request_cnt"];
                int time_limit = input_parameter["time_limit"];
                double min_fidelity = input_parameter["min_fidelity"];
                double max_fidelity = input_parameter["max_fidelity"];
                // double entangle_lambda = input_parameter["entangle_lambda"];
                // double entangle_time = input_parameter["entangle_time"];
                double entangle_prob = input_parameter["entangle_prob"];
                double swap_prob = input_parameter["swap_prob"];
                double fidelity_threshold = input_parameter["fidelity_threshold"];
                // int length_upper, length_lower;
                // if(input_parameter["path_length"] == -1) {
                //     length_upper = num_nodes;
                //     length_lower = 6;
                // } else {
                //     length_upper = input_parameter["path_length"] + 1;
                //     length_lower = input_parameter["path_length"] - 1;
                // }

                int sum_has_path = 0;
                #pragma omp parallel for
                for(int r = 0; r < round; r++) {
                    string filename = file_path + "input/round_" + to_string(r) + ".input";
                    ofstream ofs;
                    ofs.open(file_path + "log/" + path_method->get_name() + "_" + X_name + "_in_" + to_string(change_value) + "_Round_" + to_string(r) + ".log");

                    time_t now = time(0);
                    char* dt = ctime(&now);
                    cerr  << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;




                    double A = 0.25, B = 0.75, tao = input_parameter["tao"], T = 0.04, n = 2;

                    Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao, entangle_prob);

                    ofs << "--------------- in round " << r << " -------------" <<endl;
                    vector<pair<int, int>> requests;
                    for(int i = 0; i < request_cnt; i++) {
                        requests.emplace_back(default_requests[r][i]);
                    }

                    Graph path_graph = graph;
                    path_graph.increase_resources(10);
                    PathMethod *new_path_method;
                    if(path_method->get_name() == "Greedy") new_path_method = new Greedy();
                    else if(path_method->get_name() == "QCAST") new_path_method = new QCAST();
                    else if(path_method->get_name() == "REPS") new_path_method = new REPS();
                    else {
                        cerr << "unknown path method" << endl;
                        assert(false);
                    }

                    new_path_method->build_paths(path_graph, requests);
                    cout << "found path" << endl;
                    map<SDpair, vector<Path>> paths = new_path_method->get_paths();
                    map<SDpair, set<Path>> paths_st;
                    for(auto [sdpair, pathss] : paths) {
                        for(Path path : pathss) {
                            paths_st[sdpair].insert(path);
                        }
                    }

                    paths.clear();
                    for(auto [sdpair, pathss] : paths_st) {
                        for(Path path : pathss) {
                            paths[sdpair].push_back(path);
                        }
                    }

                    int path_len = 0, path_cnt = 0, mx_path_len = 0;

                    int has_path = 0;
                    for(SDpair sdpair : requests) {
                        int mi_path_len = INF;
                        has_path += !paths[sdpair].empty();
                        for(Path path : paths[sdpair]) {
                            mi_path_len = min(mi_path_len, (int)path.size());
                            for(int i = 1; i < (int)path.size(); i++) {
                                assert(graph.adj_set[path[i]].count(path[i - 1]));
                            }
                        }
                        if(mi_path_len != INF) {
                            mx_path_len = max(mx_path_len, mi_path_len);
                            path_cnt++;
                            path_len += mi_path_len;
                        }
                    }

                    sum_has_path += has_path;
                    cerr << "Path method: " << path_method->get_name() << "\n";
                    cerr << "Request cnt: " << request_cnt << "\n";
                    cerr << "Has Path cnt: " << has_path << "\n";
                    cerr << "Avg path length = " << path_len / (double)path_cnt << "\n";
                    cerr << "Max path length = " << mx_path_len << "\n";
                    vector<AlgorithmBase*> algorithms;
                    algorithms.emplace_back(new WernerAlgo(graph,requests,paths));
                    algorithms.emplace_back(new MyAlgo1(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo2(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo3(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo4(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo5(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo6(graph, requests, paths));


                    #pragma omp parallel for
                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        algorithms[i]->run();
                    }



                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        for(string Y_name : Y_names) {
                            result[r][algorithms[i]->get_name()][Y_name] = algorithms[i]->get_res(Y_name);
                        }
                    }

                    now = time(0);
                    dt = ctime(&now);
                    cerr << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;
                    ofs.close();

                    for(auto &algo : algorithms){
                        delete algo;
                    }
                    algorithms.clear();

                }

                map<string, map<string, double>> sum_res;
                // for(string algo_name : algo_names){
                //     for(int r = 0; r < round; r++){
                //         result[r][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
                //         result[r][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
                //         result[r][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
                //         result[r][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
                //         result[r][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                //     }
                // }

                for(string Y_name : Y_names) {
                    string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                    ofstream ofs;
                    ofs.open(file_path + filename, ios::app);
                    ofs << change_value << ' ';

                    for(string algo_name : algo_names){
                        for(int r = 0; r < round; r++){
                            sum_res[algo_name][Y_name] += result[r][algo_name][Y_name];
                        }
                        ofs << sum_res[algo_name][Y_name] / round << ' ';
                    }
                    ofs << endl;
                    ofs.close();
                }
            }
        }
    }
    return 0;
}