//
// Created by 罗程阳 on 2024/3/20.
//

#include <iostream>
#include <algorithm>
#include <queue>
#include <map>
#include <time.h>
#include "utilities/multilayer_graph.h"
#include "utilities/projected_graph.h"
#include "utilities/search_space.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
//string dataset_name = "stack";
string dataset_name = "higgs";
//string dataset_name = "dblp";
//string dataset_name = "homo";
//string dataset_name = "sacchcere";

int k;
int d;
int number_of_layers;
MultilayerGraph multilayer_graph;
ProjectedGraph projected_graph;

vector<int> result;

vector<int> peeling_projected_graph() {
    vector<int> bin;  // store vertices
    bin.reserve(projected_graph.maximum_vertex + 1);
    vector<int> deg(projected_graph.maximum_vertex + 1);  // the degree of each vertex
    vector<int> pos(projected_graph.maximum_vertex + 1);  // the position of each vertex in bin
    vector<int> bou(projected_graph.maximum_degree + 1);  // the boundary of each degree value in bin

    //// initialize arrays
    vector<vector<int>> degree_set(projected_graph.maximum_degree + 1);
    for (auto &vertex : projected_graph.vertex_iterator) {
        int temp_deg = projected_graph.adjacency_list[vertex].size();
        degree_set[temp_deg].push_back(vertex);
        deg[vertex] = temp_deg;
    }
    int pos_count = 0;
    for (int i = 0; i < projected_graph.maximum_degree + 1; ++i) {
        bou[i] = pos_count;
        for (auto &vertex : degree_set[i]) {
            bin.push_back(vertex);
            pos[vertex] = pos_count;
            pos_count++;
        }
    }

    //// start peeling
    for (int i = 0; i < projected_graph.maximum_vertex + 1; i++) {
        int v = bin[i];
        if (i >= bou[k + 1]) {
            break;
        }
        for (auto &u : projected_graph.adjacency_list[v]) {
            if (deg[u] > deg[v]) {
                int deg_u = deg[u];
                int pos_u = pos[u];
                int pos_w = bou[deg_u];  // w is the first vertex in the group that u is in;
                int w = bin[pos_w];
                if (u != w) {
                    // swap u and w;
                    bin[pos_u] = w;
                    bin[pos_w] = u;
                    pos[u] = pos_w;
                    pos[w] = pos_u;
                }
                bou[deg_u]++;
                deg[u]--;
            }
        }
    }

    //// get projected (k+1)-core
    vector<int> projected_result;
    for (int i = bou[k+1]; i < bin.size(); i++) {
        projected_result.push_back(bin[i]);
    }
    cout << "num of (k+1)-core:";
    cout << projected_result.size() << endl;
    return projected_result;
}

struct cmp{

    bool operator ()(const SearchNode& a, const SearchNode& b)
    {
        return a.candidate_vertex.size() < b.candidate_vertex.size();
    }
};

void prune1_delete_vertex(SearchNode &node) {
    node.get_top_d_degree_delete();
    map<int, int> deg;
    map<int, vector<int>> degree_set;
    int max_min_deg = -1;
    for (auto &vertex : node.candidate_vertex) {
        int min_deg = min(node.top_d_degree[vertex], (int)node.adjacency_list_projected[vertex].size() - 1);
        if (min_deg < 0)
            min_deg = 0;
//        int min_deg = node.top_d_degree[vertex];
        deg[vertex] = min_deg;
        degree_set[min_deg].push_back(vertex);
        max_min_deg = max(max_min_deg, min_deg);
    }

    vector<int> bin;
    map<int, int> pos;
    vector<int> bou(max_min_deg + 1);
    int pos_count = 0;
    for (int i = 0; i < max_min_deg + 1; i++) {
        bou[i] = pos_count;
        if (degree_set.find(i) != degree_set.end()) {
            for (auto vertex : degree_set[i]) {
                bin.push_back(vertex);
                pos[vertex] = pos_count;
                pos_count++;
            }
        }
    }
//    set<int> delete_num;
    for (int i = 0; i < bin.size(); i++) {
        int vertex = bin[i];
//        if (i < 100)
//            cout << vertex << endl;
        if (deg[vertex] < k) {
            set<int> update_candidate;
            for (int layer = 0; layer < number_of_layers; layer++) {
                for (auto neighbor: node.adjacency_list[vertex][layer]) {
                    if (node.adjacency_list[neighbor][layer].size() == deg[neighbor]) {
                        update_candidate.insert(neighbor);
                    }
                    if (node.adjacency_list_projected[neighbor].size() - 1 == deg[neighbor]) {
                        update_candidate.insert(neighbor);
                    }
                }
            }
            node.delete_vertex(vertex);
//            delete_num.insert(vertex);
            bou[deg[vertex]]++;
            if (deg[vertex] > 0) {
                bou[deg[vertex] - 1] = bou[deg[vertex]];
            }
            for (auto u : update_candidate) {
                int new_degree = node.get_com_degree(u);
                if (deg[u] - new_degree > 1) {
                    cout << deg[u] << " ";
                    cout << new_degree << endl;
                    cout << deg[u] - new_degree << "!!!" << endl;
                }
//                cout << deg[u] - new_degree << endl;
                if (new_degree == deg[u] - 1) {
                    int deg_u = deg[u];
                    int pos_u = pos[u];
                    int pos_w = bou[deg_u];  // w is the first vertex in the group that u is in;
                    int w = bin[pos_w];
                    if (u != w) {
                        // swap u and w;
                        bin[pos_u] = w;
                        bin[pos_w] = u;
                        pos[u] = pos_w;
                        pos[w] = pos_u;
                    }
                    bou[deg_u]++;
                    deg[u] = new_degree;
                }
            }
        } else {
            break;
        }
    }
    vector<int> new_candidate_vertex;
    if (bou.size() >= k + 1) {
        int start = bou[k];
        for (; start < bin.size(); start++) {
            new_candidate_vertex.push_back(bin[start]);
        }
    }
    node.candidate_vertex = new_candidate_vertex;
    cout << "Final:" << new_candidate_vertex.size() << endl;
    node.is_prune1 = true;
}
void prune1_composition(SearchNode &node, priority_queue<SearchNode, vector<SearchNode>, cmp> &T) {
    node.get_degree_composition(multilayer_graph, projected_graph);
    node.get_top_d_degree();
    unordered_map<int, int> deg;
    unordered_map<int, vector<int>> degree_set;
    int max_min_deg = -1;
    for (auto &vertex : node.candidate_vertex) {
        int min_deg = min(node.top_d_degree[vertex], (int)node.vertex_projected_degree[vertex] - 1);
        if (min_deg < 0)
            min_deg = 0;
        deg[vertex] = min_deg;
        degree_set[min_deg].push_back(vertex);
        max_min_deg = max(max_min_deg, min_deg);
    }

    vector<int> bin;
    unordered_map<int, int> pos;
    vector<int> bou(max_min_deg + 1);
    int pos_count = 0;
    for (int i = 0; i < max_min_deg + 1; i++) {
        bou[i] = pos_count;
        if (degree_set.find(i) != degree_set.end()) {
            for (auto vertex : degree_set[i]) {
                bin.push_back(vertex);
                pos[vertex] = pos_count;
                pos_count++;
            }
        }
    }
//    set<int> delete_num;
    for (int i = 0; i < bin.size(); i++) {
        int vertex = bin[i];
        if (deg[vertex] < k) {
            set<int> update_candidate;
            for (auto layer : node.layer_composition) {
                for (auto neighbor: multilayer_graph.adjacency_list[vertex][layer]) {
                    if (node.candidate_vertex_set.find(neighbor) != node.candidate_vertex_set.end()) {
                        if (node.vertex_layer_degree[neighbor][layer] == deg[neighbor]) {
                            update_candidate.insert(neighbor);
                        }
                        node.vertex_layer_degree[neighbor][layer]--;
                    }
                }
            }
            for (auto neighbor : projected_graph.adjacency_list[vertex]) {
                if (node.candidate_vertex_set.find(neighbor) != node.candidate_vertex_set.end()) {
                    if (node.vertex_projected_degree[neighbor] - 1 == deg[neighbor]) {
                        update_candidate.insert(neighbor);
                    }
                    node.vertex_projected_degree[neighbor]--;
                }
            }
//            for (auto neighbor : projected_graph.adjacency_list[vertex]) {
//                if (node.candidate_vertex_set.find(neighbor) != node.candidate_vertex_set.end() && node.vertex_projected_degree[neighbor] -1 == deg[neighbor]) {
//                    update_candidate.insert(neighbor);
//                }
//            }
            node.candidate_vertex_set.erase(vertex);
            bou[deg[vertex]]++;
            if (deg[vertex] > 0) {
                bou[deg[vertex] - 1] = bou[deg[vertex]];
            }
            for (auto u : update_candidate) {
                int new_degree = node.update_vertex_degree(u, vertex, multilayer_graph, projected_graph);
                if (deg[u] - new_degree > 1) {
                    cout << deg[u] << " ";
                    cout << new_degree << endl;
                    cout << deg[u] - new_degree << "???" << endl;
                }
//                cout << deg[u] - new_degree << endl;
                if (new_degree == deg[u] - 1) {
                    int deg_u = deg[u];
                    int pos_u = pos[u];
                    int pos_w = bou[deg_u];  // w is the first vertex in the group that u is in;
                    int w = bin[pos_w];
                    if (u != w) {
                        // swap u and w;
                        bin[pos_u] = w;
                        bin[pos_w] = u;
                        pos[u] = pos_w;
                        pos[w] = pos_u;
                    }
                    bou[deg_u]++;
                    deg[u] = new_degree;
                }
            }
        } else {
            break;
        }
    }
    cout << "Prune 1:" << node.candidate_vertex_set.size() << endl;
    node.candidate_vertex.clear();
    for (auto &vertex : node.candidate_vertex_set) {
        node.candidate_vertex.push_back(vertex);
    }
    node.is_prune1 = true;
    if (node.layer_composition.size() != d && node.candidate_vertex_set.size() != 0) {
        T.push(node);
    } else if (node.layer_composition.size() == d){
        if (node.candidate_vertex.size() > result.size()) {
            result = node.candidate_vertex;
        }
    }
}

void prune2(SearchNode &node, priority_queue<SearchNode, vector<SearchNode>, cmp> &T) {
//    node.get_degree_composition(multilayer_graph, projected_graph);
    node.cal_bitmap(k);
    //// put vertices into different groups
    vector<vector<int>> num_one_vertex(number_of_layers+1);
//    map<int, int> vertex_num_one;
    for (auto vertex : node.candidate_vertex_set) {
        int num_one = node.vertex_bitmap[vertex].count_one();
        num_one_vertex[num_one].push_back(vertex);
    }
    map<vector<unsigned int>, vector<int>> new_nodes;
    for (int i = d; i <= number_of_layers; i++) {
        for (int vertex : num_one_vertex[i]) {
            bool is_insert = false;
            for (auto &item : new_nodes) {
                if (node.vertex_bitmap[vertex].compare_other_bitmap(item.first)) {
                    item.second.push_back(vertex);
                    is_insert = true;
                }
            }
            if (!is_insert) {
//            if (i > d && )
                new_nodes[node.vertex_bitmap[vertex].value].push_back(vertex);
            }
        }
    }
//    map<vector<unsigned int>, vector<int>> new_nodes;
//    for (int i = d; i < number_of_layers; i++) {
//        for (int vertex : num_one_vertex[i]) {
//            bool is_insert = false;
//            for (auto &item : new_nodes) {
//                if (node.vertex_bitmap[vertex].compare_other_bitmap(item.first)) {
//                    item.second.push_back(vertex);
//                    is_insert = true;
//                }
//            }
//            if (!is_insert) {
//                new_nodes[node.vertex_bitmap[vertex].value].push_back(vertex);
//            }
//        }
//    }


//    cout << new_nodes.size() << endl;
    if (new_nodes.size() == 1) {
        auto item = new_nodes.begin();
        if (item->second.size() > result.size()) {
            result = item->second;
        }
    } else {
        for (auto &item: new_nodes) {
            vector<int> temp_composition;
            for (int i = 0; i < number_of_layers; i++) {
                int pos1 = i / 32;
                int pos2 = i % 32;
                int temp_value = item.first[pos1];
                if (((temp_value >> pos2) & 0x01) > 0) {
                    temp_composition.push_back(i);
//                    cout << i << " ";
                }
            }
//            cout << endl;
            SearchNode node(item.second, number_of_layers, d);
            node.layer_composition = temp_composition;
            if (node.layer_composition.size() == d) {
                node.is_search = true;
            }
            T.push(node);
        }
    }
}

void prune1_query_vertex_vec(SearchNode &node, priority_queue<SearchNode, vector<SearchNode>, cmp> &T, int query_vertex) {
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1, false);
    for (auto v : node.candidate_vertex) {
        vertex_exist[v] = true;
    }
    node.get_degree_composition_vec(multilayer_graph, projected_graph, vertex_exist);
    node.get_top_d_degree();
    unordered_map<int, int> deg;
    unordered_map<int, vector<int>> degree_set;
    int max_min_deg = -1;
    for (auto &vertex : node.candidate_vertex) {
        int min_deg = min(node.top_d_degree[vertex], (int)node.vertex_projected_degree[vertex] - 1);
        if (min_deg < 0)
            min_deg = 0;
        deg[vertex] = min_deg;
        degree_set[min_deg].push_back(vertex);
        max_min_deg = max(max_min_deg, min_deg);
    }

    vector<int> bin;
    unordered_map<int, int> pos;
    vector<int> bou(max_min_deg + 1);
    int pos_count = 0;
    for (int i = 0; i < max_min_deg + 1; i++) {
        bou[i] = pos_count;
        if (degree_set.find(i) != degree_set.end()) {
            for (auto vertex : degree_set[i]) {
                bin.push_back(vertex);
                pos[vertex] = pos_count;
                pos_count++;
            }
        }
    }
//    set<int> delete_num;
    for (int i = 0; i < bin.size(); i++) {
        int vertex = bin[i];
        if (deg[vertex] < k) {
            set<int> update_candidate;
            for (auto layer : node.layer_composition) {
                for (auto neighbor: multilayer_graph.adjacency_list[vertex][layer]) {
                    if (vertex_exist[neighbor]) {
                        if (node.vertex_layer_degree[neighbor][layer] == deg[neighbor]) {
                            update_candidate.insert(neighbor);
                        }
                        node.vertex_layer_degree[neighbor][layer]--;
                    }
                }
            }
            for (auto neighbor : projected_graph.adjacency_list[vertex]) {
                if (vertex_exist[neighbor]) {
                    if (node.vertex_projected_degree[neighbor] - 1 == deg[neighbor]) {
                        update_candidate.insert(neighbor);
                    }
                    node.vertex_projected_degree[neighbor]--;
                }
            }
            vertex_exist[vertex] = false;
            bou[deg[vertex]]++;
            if (deg[vertex] > 0) {
                bou[deg[vertex] - 1] = bou[deg[vertex]];
            }
            for (auto u : update_candidate) {
                int new_degree = node.update_vertex_degree(u, vertex, multilayer_graph, projected_graph);
                if (deg[u] - new_degree > 1) {
                    cout << deg[u] << " ";
                    cout << new_degree << endl;
                    cout << deg[u] - new_degree << "???" << endl;
                }
//                cout << deg[u] - new_degree << endl;
                if (new_degree == deg[u] - 1) {
                    int deg_u = deg[u];
                    int pos_u = pos[u];
                    int pos_w = bou[deg_u];  // w is the first vertex in the group that u is in;
                    int w = bin[pos_w];
                    if (u != w) {
                        // swap u and w;
                        bin[pos_u] = w;
                        bin[pos_w] = u;
                        pos[u] = pos_w;
                        pos[w] = pos_u;
                    }
                    bou[deg_u]++;
                    deg[u] = new_degree;
                }
            }
        } else {
            break;
        }
    }
    if (vertex_exist[query_vertex]) {
        vector<int> new_candidate_vertex;
        for (auto vertex : node.candidate_vertex) {
            if (vertex_exist[vertex]) {
                new_candidate_vertex.push_back(vertex);
            }
        }
        node.candidate_vertex = new_candidate_vertex;
        cout << "Prune 1:" << node.candidate_vertex.size() << endl;
        node.is_prune1 = true;
        if (node.layer_composition.size() != d && node.candidate_vertex.size() != 0) {
            T.push(node);
        } else if (node.layer_composition.size() == d) {
            vector<int> temp_result;
            queue<int> Q;
            Q.push(query_vertex);
            vector<bool> vertex_visited(multilayer_graph.maximum_vertex + 1, false);
            while (!Q.empty()) {
                int front_vertex = Q.front();
                Q.pop();
                temp_result.push_back(front_vertex);
                for (auto neighbor : projected_graph.adjacency_list[front_vertex]) {
                    if (vertex_visited[neighbor])
                        continue;
                    vertex_visited[neighbor] = true;
                    if (vertex_exist[neighbor]) {
                        Q.push(neighbor);
                    }
                }
            }
            if (temp_result.size() > result.size()) {
                result = temp_result;
            }
        }
    } else {
        cout << "This node does not contain the query vertex!" << endl;
    }
}

void prune2_query_vertex_vec(SearchNode &node, priority_queue<SearchNode, vector<SearchNode>, cmp> &T, int query_vertex) {
//    node.get_degree_composition(multilayer_graph, projected_graph);
    node.cal_bitmap(k);
    //// put vertices into different groups
    vector<vector<int>> num_one_vertex(number_of_layers+1);
//    map<int, int> vertex_num_one;
    for (auto vertex : node.candidate_vertex) {
        int num_one = node.vertex_bitmap[vertex].count_one();
        num_one_vertex[num_one].push_back(vertex);
    }
    map<vector<unsigned int>, vector<int>> new_nodes;
    for (int i = d; i <= number_of_layers; i++) {
        for (int vertex : num_one_vertex[i]) {
            bool is_insert = false;
            for (auto &item : new_nodes) {
                if (node.vertex_bitmap[vertex].compare_other_bitmap(item.first)) {
                    item.second.push_back(vertex);
                    is_insert = true;
                }
            }
            if (!is_insert) {
                new_nodes[node.vertex_bitmap[vertex].value].push_back(vertex);
            }
        }
    }
//    map<vector<unsigned int>, vector<int>> new_nodes;
//    for (int i = d; i < number_of_layers; i++) {
//        for (int vertex : num_one_vertex[i]) {
//            bool is_insert = false;
//            for (auto &item : new_nodes) {
//                if (node.vertex_bitmap[vertex].compare_other_bitmap(item.first)) {
//                    item.second.push_back(vertex);
//                    is_insert = true;
//                }
//            }
//            if (!is_insert) {
//                new_nodes[node.vertex_bitmap[vertex].value].push_back(vertex);
//            }
//        }
//    }


//    cout << new_nodes.size() << endl;
    if (new_nodes.size() == 1) {
        auto item = new_nodes.begin();
        vector<int> candidate_vertex;
        vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1, false);
        for (auto v : item->second) {
            vertex_exist[v] = true;
        }
        if (!vertex_exist[query_vertex]) {
            return;
        }
        vector<int> temp_result;
        queue<int> Q;
        Q.push(query_vertex);
        vector<bool> vertex_visited(multilayer_graph.maximum_vertex + 1, false);
        while (!Q.empty()) {
            int front_vertex = Q.front();
            Q.pop();
            temp_result.push_back(front_vertex);
            for (auto neighbor : projected_graph.adjacency_list[front_vertex]) {
                if (vertex_visited[neighbor])
                    continue;
                vertex_visited[neighbor] = true;
                if (vertex_exist[neighbor]) {
                    Q.push(neighbor);
                }
            }
        }
//        while (!Q.empty()) {
//            int front_vertex = Q.front();
//            Q.pop();
//            temp_result.push_back(front_vertex);
//            for (auto neighbor : projected_graph.adjacency_list[front_vertex]) {
//                if (vertex_visited.find(neighbor) != vertex_visited.end())
//                    continue;
//                vertex_visited.insert(neighbor);
//                if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
//                    Q.push(neighbor);
//                }
//            }
//        }
        if (temp_result.size() > result.size()) {
            result = temp_result;
        }
    } else {
        for (auto &item: new_nodes) {
            if (std::find(item.second.begin(), item.second.end(),query_vertex) == item.second.end()) {
                continue;
            }
            vector<int> temp_composition;
            for (int i = 0; i < number_of_layers; i++) {
                int pos1 = i / 32;
                int pos2 = i % 32;
                int temp_value = item.first[pos1];
                if (((temp_value >> pos2) & 0x01) > 0) {
                    temp_composition.push_back(i);
//                    cout << i << " ";
                }
            }
//            cout << endl;
            SearchNode node;
            node.initialize_vec(item.second, number_of_layers, d);
            node.layer_composition = temp_composition;
            if (node.layer_composition.size() == d) {
                node.is_search = true;
            }
            T.push(node);
        }
    }
}

void MLCS_query_vertex_vec(int query_vertex) {
    //// search connected (k+1)-core in the projected graph
    vector<int>projected_k1_core = peeling_projected_graph();
//    vector<int>projected_k1_core(multilayer_graph.maximum_vertex + 1);
//    for (int i = 0; i < projected_k1_core.size(); i++) {
//        projected_k1_core[i] = i;
//    }
    priority_queue<SearchNode, vector<SearchNode>, cmp> T;
//    SearchNode root(projected_k1_core, number_of_layers, d);
    SearchNode root;
    root.initialize_vec(projected_k1_core, number_of_layers, d);
//    root.get_degree(multilayer_graph, projected_graph);
    root.layer_composition.resize(number_of_layers);
    for (int i = 0; i < number_of_layers; i++) {
        root.layer_composition[i] = i;
    }
    T.push(root);
    while (!T.empty()) {
        auto node = T.top();
        if (node.candidate_vertex.size() < result.size())
            break;
        T.pop();
        if (!node.is_prune1) {
            prune1_query_vertex_vec(node, T, query_vertex);
        } else {
            prune2_query_vertex_vec(node, T, query_vertex);
        }
    }
    cout << "Final result num:" << result.size() << endl;
}


void exp() {
    ifstream test_case_file("datasets/" + dataset_name + "/" + dataset_name + "_test_case.txt");
    ofstream exp_file("datasets/" + dataset_name + "/" + dataset_name + "_prune_exp.txt");
    string line;
    while (getline(test_case_file, line)) {
        if (line.size() == 0)
            continue;
        istringstream iss(line);
        vector<string> words;
        string word;
        while (iss >> word) {
            words.push_back(word);
        }
        if (word.size() == 0)
            continue;
        result.clear();
        k = stoi(words[0]);
        d = stoi(words[1]);
        int query_vertex = stoi(words[2]);
        exp_file << "k:" << k << " " << "d:" << d << endl;
        cout << "Prune on " << dataset_name << endl;
        cout << "k:" << k << " " << "d:" << d << endl;
        cout << "Query vertex :" << query_vertex << endl;
        exp_file << "Query vertex :" << query_vertex << endl;
        auto start = std::chrono::high_resolution_clock::now();
//        MLCS_query_vertex(query_vertex);
        MLCS_query_vertex_vec(query_vertex);
//        MLCS_query_vertex_core(query_vertex);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << "Result Num:" << result.size() << endl;
        exp_file << "Result Num:" << result.size() << endl;
        cout << "Search Time (s): " << duration.count() / 1000.0 << std::endl;
        exp_file << "Search Time (s): " << duration.count() / 1000.0 <<  endl << endl;
    }
    exp_file.close();
}

int main(int argc, char* argv[]) {
    srand(time(0));
    clock_t start_time, end_time;
    dataset_name = argv[1];
    string dataset_path = "datasets/" + dataset_name + "/" + dataset_name + ".txt";
    cout << "Loading graph...." << endl;
    // load dataset
    MultilayerGraph multilayer_graph_dataset(dataset_path);
    multilayer_graph = multilayer_graph_dataset;
    ProjectedGraph projected_graph_dataset(multilayer_graph);
    projected_graph = projected_graph_dataset;
    number_of_layers = multilayer_graph.number_of_layers;

    exp();

    return 0;
}
