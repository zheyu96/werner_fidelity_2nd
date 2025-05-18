#ifndef __NODE_H
#define __NODE_H

#include "../../config.h"

using namespace std;

class Node {
    int id;
    int memory, time_limit;
    double swap_prob;
    vector<int> memory_table;
    vector<int> neighbor;
public:
    Node(int id, int memory, int time_limit, double swap_prob);
    ~Node();
    int get_id();
    int get_memory();
    int get_memory_at(int t);
    double get_swap_prob();
    vector<int> get_neighbot();
    void add_neighbor(int neighbor_id);
    void reserve_memory(int t, int amount);
    void reserve_memory(int t);
};

#endif