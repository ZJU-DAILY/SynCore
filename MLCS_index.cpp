//
// Created by 罗程阳 on 2024/4/8.
//

#include <iostream>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <time.h>
#include <chrono>
#include <unordered_set>
#include "utilities/multilayer_graph.h"
#include "utilities/projected_graph.h"
#include "utilities/search_space.h"
#include "utilities/trie.h"
#include <fstream>
#include <sstream>
#include <bitset>
#include <stack>

using namespace std;

//string dataset_name = "higgs";
//string dataset_name = "dblp";
//string dataset_name = "homo";
//string dataset_name = "slashdot";
//string dataset_name = "obamainisrael";
string dataset_name = "yeast";
//string dataset_name = "sacchcere";
//string dataset_name = "stack";
//string dataset_name = "wiki";

int number_of_layers;
MultilayerGraph multilayer_graph;
ProjectedGraph projected_graph;
int k, d;
int k_max = 0;
vector<vector<vector<vector<vector<int>>>>> PMC_index;
vector<vector<vector<vector<BitMap>>>> PMC_index_bitmap;
vector<int> result;
vector<trie_node*> tries;
vector<vector<vector<trie_node*>>> trie_vertex_node;

double graph_size = 0;
double index_size = 0;
double build_time = 0;

void load_index() {
    PMC_index.resize(multilayer_graph.maximum_vertex + 1);
    ifstream index_file("datasets/" + dataset_name + "/" + dataset_name + "_index.txt");
    string line;
    int current_v = -1;
    int current_k = -1;
    while (getline(index_file, line)) {
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
        if (words[0] == "v") {
            current_v = stoi(words[1]);
            PMC_index[current_v].push_back({});
        } else if ((words[0] == "k")) {
            current_k = stoi(words[1]);
            k_max = max(k_max, current_k);
            PMC_index[current_v].push_back({});
        } else {
            vector<int> layers;
            for (auto str : words) {
                layers.push_back(stoi(str));
            }
            if (PMC_index[current_v][current_k].size() == 0) {
                PMC_index[current_v][current_k].resize(number_of_layers + 1);
            }
            PMC_index[current_v][current_k][layers.size()].push_back(layers);
        }
    }
    index_file.close();
}

void convert_bitmap_index() {
    long long index_size = 0;
    long long graph_size = 0;
    PMC_index_bitmap.resize(multilayer_graph.maximum_vertex + 1);
    index_size += (multilayer_graph.maximum_vertex + 1) * 4;
    graph_size += (multilayer_graph.maximum_vertex + 1) * number_of_layers * 4;
    long long bitmap_num = 0;

    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        PMC_index_bitmap[v].resize(PMC_index[v].size());
        index_size += PMC_index[v].size() * 4;
        for (int temp_k = 1; temp_k < PMC_index[v].size(); temp_k++) {
            PMC_index_bitmap[v][temp_k].resize(PMC_index[v][temp_k].size());
            for (int temp_d = 1; temp_d < PMC_index[v][temp_k].size(); temp_d++) {
                for (auto layers : PMC_index[v][temp_k][temp_d]) {
                    bitmap_num++;
                    BitMap temp_bitmap(number_of_layers);
                    for (auto layer : layers) {
                        temp_bitmap.set_one(layer);
                    }
                    PMC_index_bitmap[v][temp_k][temp_d].push_back(temp_bitmap);
                }
            }
        }
    }
    index_size += (bitmap_num * number_of_layers) / 8;
    cout << "Index Size:" << (double)index_size / 1024.0 / 1024.0 << "MB" << endl;

    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        for (auto l = 0; l < number_of_layers; l++) {
            graph_size += multilayer_graph.adjacency_list[v][l].size() * 4;
        }
    }
    cout << "Graph Size:" << (double)graph_size / 1024.0 / 1024.0 << "MB" << endl;


//    //// write index using bitmap
//    ofstream index_file_bitmap("datasets/" + dataset_name + "/" + dataset_name + "_index_bitmap", std::ios::binary | std::ios::out);
//    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
//        cout << "v " << v << endl;
//        char ch = 'v';
//        index_file_bitmap.write(&ch, sizeof(char));
//        index_file_bitmap.write(reinterpret_cast<const char*>(&v), sizeof(int));
//        for (int k = 1; k < PMC_index_bitmap[v].size(); k++) {
//            cout << "k " << k << endl;
//            ch = 'k';
//            index_file_bitmap.write(&ch, sizeof(char));
//            index_file_bitmap.write(reinterpret_cast<const char*>(&k), sizeof(int));
//            for (int d = 1; d < PMC_index_bitmap[v][k].size(); d++) {
//                if (PMC_index_bitmap[v][k][d].size() == 0)
//                    continue;
////                cout << "d=" << d << endl;
////                index_file << "d=" << d << endl;
//                for (auto bitmap : PMC_index_bitmap[v][k][d]) {
//                    index_file_bitmap.write(reinterpret_cast<const char*>(bitmap.value.data()), bitmap.value.size()*sizeof(int));
//                }
//            }
//        }
//    }
//
//    index_file_bitmap.close();

}

void index_correct_test() {
//    //// construct the projected graph of the multilayer graph
//    ProjectedGraph projected_graph(multilayer_graph);

//    vector<vector<int>> layer_composition = {{0}, {1}, {3,4}, {2,7}, {1, 5, 9}};
    vector<vector<int>> layer_composition = {{0},
                                             {1},
                                             {2},
                                             {3},
                                             {1, 2},
                                             {0, 3},
                                             {1, 2, 3},
                                             {0, 2},
                                             {0, 2, 3}};
//    vector<vector<int>> layer_composition = {{0}, {1}, {5}, {7}, {3,4}, {2,7}};

    int k = 1;

    int final_num = 0;
    for (auto &composition: layer_composition) {
        for (auto &num: composition) {
            cout << num << " ";
        }
        cout << endl;

        // test correctness
        bool change = true;
        auto test_multilayer_graph = multilayer_graph;
        auto test_projected_graph = projected_graph;
        vector<bool> test_exist(multilayer_graph.maximum_vertex + 1, true);
//        int left_count = multilayer_graph.maximum_vertex + 1;
        int left_count = multilayer_graph.maximum_vertex + 1;
        while (change) {
            change = false;
            for (int i = 0; i < multilayer_graph.maximum_vertex + 1; i++) {
                bool is_delete = false;
                if (test_exist[i]) {
                    int n1 = 0;
                    for (auto neighbor: test_projected_graph.adjacency_list[i]) {
                        if (test_exist[neighbor]) {
                            ++n1;
                        }
                    }
                    if (n1 < k + 1) {
                        is_delete = true;
                    }
                    for (auto layer: composition) {
                        int n2 = 0;
                        for (auto neighbor: test_multilayer_graph.adjacency_list[i][layer]) {
                            if (test_exist[neighbor]) {
                                ++n2;
                            }
                        }
                        if (n2 < k)
                            is_delete = true;
                    }

                    if (is_delete) {
                        change = true;
                        test_exist[i] = false;
                        left_count--;
                    }
                }
            }
            if (!change)
                break;
        }
        cout << "correct count:" << left_count << endl;

//        int count_num = 0;
//        for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
//            if (vertex_core_number[v].find(composition) != vertex_core_number[v].end()) {
//                if (vertex_core_number[v][composition] >= k) {
//                    count_num++;
//                }
//            }
//        }
//        cout << "decomposition count:" << count_num << endl;
        int count_num = 0;
        for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
            if (k > PMC_index[v].size() - 1 || PMC_index[v][k].size() == 0) {
                continue;
            }
            bool exist = false;
            for (int layer_num = composition.size(); layer_num < number_of_layers + 1; layer_num++) {
                if (exist || layer_num > PMC_index[v][k].size() - 1)
                    break;
                for (auto layers : PMC_index[v][k][layer_num]) {
                    bool dominate = true;
                    for (auto layer : composition) {
                        if (std::find(layers.begin(), layers.end(), layer) == layers.end()) {
                            dominate = false;
                            break;
                        }
                    }
                    if (dominate) {
                        exist = true;
                        break;
                    }
                }
            }
            if (exist) {
                count_num++;
            }
        }
        cout << "index count:" << count_num << endl;
    }
}

void generate_layer_composition(int layer, int select, vector<vector<int>> &layer_composition, vector<int> &layer_array, vector<bool> &layer_visited) {
    if (select == d) {
        vector<int> temp_layer_composition;
        for (int i = 0; i < number_of_layers; ++i) {
            if (layer_visited[i]) {
                temp_layer_composition.push_back(i);
            }
        }
        layer_composition.push_back(temp_layer_composition);
        return;
    }
    if (layer == number_of_layers)
        return;
    layer_visited[layer] = true;
    generate_layer_composition(layer + 1, select + 1, layer_composition, layer_array, layer_visited);
    layer_visited[layer] = false;
    generate_layer_composition(layer + 1, select, layer_composition, layer_array, layer_visited);
}

void generate_sublayer_composition(int layer, int select, vector<int> &layers, vector<vector<int>> &sub_layers, vector<bool> &layer_visited) {
    if (select == d) {
        vector<int> temp_layer_composition;
        for (int i = 0; i < number_of_layers; ++i) {
            if (layer_visited[i]) {
                temp_layer_composition.push_back(i);
            }
        }
        sub_layers.push_back(temp_layer_composition);
        return;
    }
    if (layer == layers.size())
        return;
    layer_visited[layers[layer]] = true;
    generate_sublayer_composition(layer + 1, select + 1, layers, sub_layers, layer_visited);
    layer_visited[layers[layer]] = false;
    generate_sublayer_composition(layer + 1, select, layers, sub_layers, layer_visited);
}

vector<vector<int>> generate_sub_layers(vector<int> &layers) {
    vector<vector<int>> sub_layers;
    if (layers.size() == d) {
        sub_layers.push_back(layers);
        return sub_layers;
    }
    vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
    generate_sublayer_composition(0, 0, layers, sub_layers, layer_visited);

    return sub_layers;
}

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

    // test correctness
    bool change = true;
    auto test_projected_graph = projected_graph;
    vector<bool> test_exist(projected_graph.maximum_vertex + 1, true);
    int left_count = projected_graph.maximum_vertex + 1;
    while (change) {
        change = false;
        for (int i = 0; i < projected_graph.maximum_vertex + 1; i++) {
            bool is_delete = false;
            if (test_exist[i]) {
                if (test_projected_graph.adjacency_list[i].size() < k + 1) {
                    is_delete = true;
                }
                if (is_delete) {
                    change = true;
                    test_projected_graph.delete_vertex(i);
                    test_exist[i] = false;
                    left_count--;
                }
            }
        }
        if (!change)
            break;
    }

    cout << "left:" << left_count << endl;

    return projected_result;
}

vector<int> baseline() {
    vector<int> result;
    //// search connected (k+1)-core in the projected graph
    vector<int>projected_k1_core = peeling_projected_graph();

    //// get the composition of d layers
    vector<vector<int>> layer_composition;
    vector<int> layer_array = multilayer_graph.layers_iterator;
    vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
    generate_layer_composition(0, 0, layer_composition, layer_array, layer_visited);

    unordered_set<int> candidate_vertex;
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1, false);
    for (auto &vertex : projected_k1_core) {
        candidate_vertex.insert(vertex);
        vertex_exist[vertex] = true;
    }
    for (auto &composition : layer_composition) {
        for (auto &num : composition) {
            cout << num << " ";
        }
        cout << endl;
        auto temp_candidate_vertex = candidate_vertex;
        auto temp_vertex_exist = vertex_exist;

        // compute minimum degree of vertex in L layers
        map<int, vector<int>> degree_set;
        vector<int> deg(multilayer_graph.maximum_vertex + 1);
        int max_min_deg = 0;
//        unordered_map<int, vector<int>>
        for (auto &vertex : candidate_vertex) {
            int min_degree = multilayer_graph.maximum_vertex - 1;
            for (auto layer : composition) {
                int temp_layer_deg = 0;
                for (auto neighbor : multilayer_graph.adjacency_list[vertex][layer]) {
//                    if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                        temp_layer_deg++;
//                    }
                    if (vertex_exist[neighbor]) {
                        temp_layer_deg++;
                    }
                }
                min_degree = min(min_degree, temp_layer_deg);
            }
            int temp_projected_degree = 0;
            for (auto neighbor : projected_graph.adjacency_list[vertex]) {
//                if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                    temp_projected_degree++;
//                }
                if (vertex_exist[neighbor]) {
                    temp_projected_degree++;
                }
            }
            min_degree = min(min_degree, temp_projected_degree - 1);
            if (min_degree < 0)
                min_degree = 0;
            deg[vertex] = min_degree;
            degree_set[min_degree].push_back(vertex);
            max_min_deg = max(max_min_deg, min_degree);
        }

        vector<int> bin, pos, bou;
        pos.resize(multilayer_graph.maximum_vertex + 1);  // the index of v in the bin
        bou.resize(max_min_deg + 1);  // the start position of bou structure
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

        // start peeling
        int i = 0;
        for (i = 0; i < bin.size(); i++) {
            if (temp_candidate_vertex.size() < result.size())
                break;
//            if (i % 1000 == 0) {
//                cout << i << "/" << bin.size() << endl;
//            }
            int vertex = bin[i];
            if (deg[vertex] < k) {
                unordered_set<int> affect_neighbor;
//                for (auto layer : composition) {
//                    for (auto neighbor : multilayer_graph.adjacency_list[vertex][layer]) {
//                        if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                            affect_neighbor.insert(neighbor);
//                        }
//                    }
//                }
                for (auto neighbor : projected_graph.adjacency_list[vertex]) {
//                    if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                        affect_neighbor.insert(neighbor);
//                    }
                    if (temp_vertex_exist[neighbor]) {
                        affect_neighbor.insert(neighbor);
                    }
                }

                bou[deg[vertex]]++;
                if (deg[vertex] > 0) {
                    bou[deg[vertex] - 1] = bou[deg[vertex]];
                }
                temp_candidate_vertex.erase(vertex);
                temp_vertex_exist[vertex] = false;
                for (auto u : affect_neighbor) {
                    int min_degree = multilayer_graph.number_of_vertex - 1;
                    for (auto layer : composition) {
                        int temp_layer_deg = 0;
                        for (auto neighbor : multilayer_graph.adjacency_list[u][layer]) {
//                            if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                                temp_layer_deg++;
//                            }
                            if (temp_vertex_exist[neighbor]) {
                                temp_layer_deg++;
                            }
                        }
                        min_degree = min(min_degree, temp_layer_deg);
                    }
                    int temp_projected_degree = 0;
                    for (auto neighbor : projected_graph.adjacency_list[u]) {
//                        if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                            temp_projected_degree++;
//                        }
                        if (temp_vertex_exist[neighbor]) {
                            temp_projected_degree++;
                        }
                    }
                    min_degree = min(min_degree, temp_projected_degree - 1);
                    if (min_degree < 0)
                        min_degree = 0;
//                    if (deg[u] - min_degree > 1)
//                        cout << "!!!" << endl;
                    if (min_degree != deg[u]) {
                        int w = bin[bou[deg[u]]];
                        int pos_w = pos[w];
                        int pos_u = pos[u];
                        bin[pos_u] = w;
                        bin[pos_w] = u;
                        pos[u] = pos_w;
                        pos[w] = pos_u;
                        bou[deg[u]]++;
                        deg[u] = min_degree;
                    }
                }
            } else {
                break;
            }
        }
//        cout << "(k+1)-core: " << temp_result.size() - i<< endl;
        cout << "temp result: " << temp_candidate_vertex.size() << endl;
        if (temp_candidate_vertex.size() > result.size()) {
            result.clear();
            for (i; i < bin.size(); i++) {
                result.push_back(bin[i]);
            }
        }
    }
    return result;
}

void MLCS(int query_vertex) {
    if (!(PMC_index[query_vertex].size() >= k + 1 && PMC_index[query_vertex][k].size() == (number_of_layers + 1))) {
        cout << "No result!" << endl;
        return;
    }
    set<vector<int>> candidate_layers_set;
    for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
        for (auto layers : PMC_index[query_vertex][k][layer_num]) {
            auto sub_layers = generate_sub_layers(layers);
            for (auto sub_layer : sub_layers) {
                candidate_layers_set.insert(sub_layer);
            }
        }
    }
    for (auto layers : candidate_layers_set) {
        vector<int> temp_result;
        unordered_set<int> visited_vertex;
        queue<int> Q;
        Q.push(query_vertex);
        while (!Q.empty()) {
            int front_vertex = Q.front();
            Q.pop();
            temp_result.push_back(front_vertex);
            for (auto neighbor: projected_graph.adjacency_list[front_vertex]) {
                if (visited_vertex.find(neighbor) != visited_vertex.end()) {
                    continue;
                }
                visited_vertex.insert(neighbor);
                if (!(PMC_index[neighbor].size() >= k + 1 &&
                      PMC_index[neighbor][k].size() == number_of_layers + 1)) {
                    continue;
                }
                bool dominate = false;
                for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
                    if (dominate) {
                        break;
                    }
                    for (auto dominate_layers: PMC_index[neighbor][k][layer_num]) {
                        bool temp_dominate = true;
                        for (auto num: layers) {
                            if (std::find(dominate_layers.begin(), dominate_layers.end(), num) ==
                                dominate_layers.end()) {
                                temp_dominate = false;
                                break;
                            }
                        }
                        if (temp_dominate) {
                            dominate = true;
                            break;
                        }
                    }
                }
                if (dominate) {
                    Q.push(neighbor);
                }
            }
        }
        if (temp_result.size() > result.size())
            result = temp_result;
    }
    cout << result.size() << endl;
}

void MLCS_bitmap(int query_vertex) {
    if (!(PMC_index[query_vertex].size() >= k + 1 && PMC_index[query_vertex][k].size() == (number_of_layers + 1))) {
        cout << "No result!" << endl;
        return;
    }
    set<vector<int>> candidate_layers_set;
    for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
        for (auto layers : PMC_index[query_vertex][k][layer_num]) {
            auto sub_layers = generate_sub_layers(layers);
            for (auto sub_layer : sub_layers) {
                candidate_layers_set.insert(sub_layer);
            }
        }
    }
    vector<bool> visited_vertex_vec(multilayer_graph.maximum_vertex + 1, false);
    for (auto layers : candidate_layers_set) {
        std::fill(visited_vertex_vec.begin(), visited_vertex_vec.end(), false);
        vector<int> temp_result;
        queue<int> Q;
        Q.push(query_vertex);
        BitMap target_bitmap(number_of_layers);
        for (auto layer : layers) {
            target_bitmap.set_one(layer);
        }
        visited_vertex_vec[query_vertex] = true;
        while (!Q.empty()) {
            int front_vertex = Q.front();
            Q.pop();
            temp_result.push_back(front_vertex);
            for (auto neighbor: projected_graph.adjacency_list[front_vertex]) {
                if (visited_vertex_vec[neighbor]) {
                    continue;
                }
                visited_vertex_vec[neighbor] = true;
                if (!(PMC_index_bitmap[neighbor].size() >= k + 1 &&
                      PMC_index_bitmap[neighbor][k].size() == number_of_layers + 1)) {
                    continue;
                }
                bool dominate = false;
                for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
                    if (dominate) {
                        break;
                    }
                    for (auto temp_bitmap: PMC_index_bitmap[neighbor][k][layer_num]) {
                        bool temp_dominate = true;
                        for (int i = 0; i < target_bitmap.value.size(); i++) {
                            unsigned int temp_bits = ((target_bitmap.value[i] & temp_bitmap.value[i]) ^ target_bitmap.value[i]);
                            if (temp_bits != 0) {
                                temp_dominate = false;
                                break;
                            }
                        }
                        if (temp_dominate) {
                            dominate = true;
                            break;
                        }
                    }
                }
                if (dominate) {
                    Q.push(neighbor);
                }
            }
        }
        if (temp_result.size() > result.size())
            result = temp_result;
    }

//    for (auto layers : candidate_layers_set) {
//        vector<int> temp_result;
//        unordered_set<int> visited_vertex;
//        queue<int> Q;
//        Q.push(query_vertex);
//        BitMap target_bitmap(number_of_layers);
//        for (auto layer : layers) {
//            target_bitmap.set_one(layer);
//        }
//        while (!Q.empty()) {
//            int front_vertex = Q.front();
//            Q.pop();
//            temp_result.push_back(front_vertex);
//            for (auto neighbor: projected_graph.adjacency_list[front_vertex]) {
//                if (visited_vertex.find(neighbor) != visited_vertex.end()) {
//                    continue;
//                }
//                visited_vertex.insert(neighbor);
//                if (!(PMC_index_bitmap[neighbor].size() >= k + 1 &&
//                        PMC_index_bitmap[neighbor][k].size() == number_of_layers + 1)) {
//                    continue;
//                }
//                bool dominate = false;
//                for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
//                    if (dominate) {
//                        break;
//                    }
//                    for (auto temp_bitmap: PMC_index_bitmap[neighbor][k][layer_num]) {
//                        bool temp_dominate = true;
//                        for (int i = 0; i < target_bitmap.value.size(); i++) {
//                            unsigned int temp_bits = ((target_bitmap.value[i] & temp_bitmap.value[i]) ^ target_bitmap.value[i]);
//                            if (temp_bits != 0) {
//                                temp_dominate = false;
//                                break;
//                            }
//                        }
//                        if (temp_dominate) {
//                            dominate = true;
//                            break;
//                        }
//                    }
//                }
//                if (dominate) {
//                    Q.push(neighbor);
//                }
//            }
//        }
//        if (temp_result.size() > result.size())
//            result = temp_result;
//    }
    cout << result.size() << endl;
}

void result_test(){
    vector<pair<int, int>> k_d_list = {{5,  1},
                                       {5,  2},
                                       {5,  3},
                                       {5,  4},
                                       {5,  5},
                                       {5,  6},
                                       {10, 1},
                                       {10, 2},
                                       {10, 3},
                                       {10, 4},
                                       {10, 5},
                                       {15, 1},
                                       {15, 2},
                                       {15, 3},
                                       {15, 4},
                                       {20, 1},
                                       {20, 2},
                                       {20, 3}};
    for (auto k_d_pair: k_d_list) {
        k = k_d_pair.first;
        d = k_d_pair.second;
        //// get the composition of d layers
        vector<vector<int>> layer_composition;
        vector<int> layer_array = multilayer_graph.layers_iterator;
        vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
        generate_layer_composition(0, 0, layer_composition, layer_array, layer_visited);
        int max_result = 0;
        for (auto composition: layer_composition) {
            int count_num = 0;
            vector<int> temp_result;
            for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                if (k > PMC_index[v].size() - 1 || PMC_index[v][k].size() == 0) {
                    continue;
                }
                bool exist = false;
                for (int layer_num = composition.size(); layer_num < number_of_layers + 1; layer_num++) {
                    if (exist || layer_num > PMC_index[v][k].size() - 1)
                        break;
                    for (auto layers: PMC_index[v][k][layer_num]) {
                        bool dominate = true;
                        for (auto layer: composition) {
                            if (std::find(layers.begin(), layers.end(), layer) == layers.end()) {
                                dominate = false;
                                break;
                            }
                        }
                        if (dominate) {
                            exist = true;
                            break;
                        }
                    }
                }
                if (exist) {
                    count_num++;
                    temp_result.push_back(v);
                }
            }
            if (count_num > max_result) {
                max_result = count_num;
                result = temp_result;
            }
//                cout << "index count:" << count_num << endl;
        }
//        cout << "k:" << k << " " << "d:" << d << endl;
//        cout << "max num:" << max_result << endl;
//        cout << result.size() << endl;
//        cout << "Query result:" << result[0] << endl;
        cout << k << " " << d << " " << result[0] << endl;
    }
//    vector<int> k_list = {5, 10, 15, 20};
//    vector<int> d_list = {1, 2, 3, 4, 5, 6, 7};
//    for (int i = 0; i < k_list.size(); i++) {
//        for (int j = 0; j < d_list.size(); j++) {
//            k = k_list[i];
//            d = d_list[j];
//            //// get the composition of d layers
//            vector<vector<int>> layer_composition;
//            vector<int> layer_array = multilayer_graph.layers_iterator;
//            vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
//            generate_layer_composition(0, 0, layer_composition, layer_array, layer_visited);
//            int max_result = 0;
//            for (auto composition : layer_composition) {
//                int count_num = 0;
//                for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
//                    if (k > PMC_index[v].size() - 1 || PMC_index[v][k].size() == 0) {
//                        continue;
//                    }
//                    bool exist = false;
//                    for (int layer_num = composition.size(); layer_num < number_of_layers + 1; layer_num++) {
//                        if (exist || layer_num > PMC_index[v][k].size() - 1)
//                            break;
//                        for (auto layers: PMC_index[v][k][layer_num]) {
//                            bool dominate = true;
//                            for (auto layer: composition) {
//                                if (std::find(layers.begin(), layers.end(), layer) == layers.end()) {
//                                    dominate = false;
//                                    break;
//                                }
//                            }
//                            if (dominate) {
//                                exist = true;
//                                break;
//                            }
//                        }
//                    }
//                    if (exist) {
//                        count_num++;
//                    }
//                }
//                if (count_num > max_result) {
//                    max_result = count_num;
//                }
////                cout << "index count:" << count_num << endl;
//            }
//            cout << "k:" << k << " " << "d:" << d << endl;
//            cout << "max num:" << max_result << endl;
//        }
//    }
}

void build_trie() {
    tries.resize(k_max + 1);
    trie_vertex_node.resize(multilayer_graph.maximum_vertex + 1);
    for (int i = 0; i < multilayer_graph.maximum_vertex + 1; i++) {
        trie_vertex_node[i].resize(k_max + 1);
    }
    long long total_size_trie = 0;
    int node_num = 0;
    int vertex_num = 0;
    for (int temp_k = 1; temp_k <= k_max; temp_k++) {
        auto *root = new trie_node(-1);
        node_num++;
        for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
            if (PMC_index[v].size() < temp_k + 1)
                continue;
            for (int temp_d = 1; temp_d < PMC_index[v][temp_k].size(); temp_d++) {
                for (auto layers : PMC_index[v][temp_k][temp_d]) {
                    trie_node* node = root;
                    for (int i = 0; i < layers.size(); i++) {
                        int layer = layers[i];
                        bool exist = false;
                        for (auto child : node->child_node) {
                            if (child->layer_id == layer) {
                                node = child;
                                exist = true;
                                break;
                            }
                        }
                        if (!exist) {
                            node_num++;
                            total_size_trie += 4;
                            total_size_trie += 4;
                            total_size_trie += 1;
                            trie_node *new_node = new trie_node(layer);
                            new_node->layer_id = layer;
                            int insert_pos = 0;
                            while (insert_pos < node->child_node.size()) {
                                if (layer > node->child_node[insert_pos]->layer_id) {
                                    insert_pos++;
                                } else {
                                    break;
                                }
                            }
                            node->child_node.insert(node->child_node.begin() + insert_pos, new_node);
                            total_size_trie += 8;
                            total_size_trie += 8;
                            new_node->father_node = node;
                            node = new_node;
                            // update the height
                            trie_node* temp_node = node;
                            while (temp_node != root) {
                                if (temp_node->father_node->height < temp_node->height + 1) {
                                    temp_node->father_node->height = temp_node->height + 1;
                                    temp_node = temp_node->father_node;
                                } else {
                                    break;
                                }
                            }
                        }
                        if (i == layers.size() - 1) {
                            trie_vertex_node[v][temp_k].push_back(node);
                            total_size_trie += 4;
                            node->vertex_vec.push_back(v);
                            vertex_num++;
                            node->dominate = true;
                        }
                    }
                }
            }
        }
        tries[temp_k] = root;
    }
    cout << "Number of Trie Nodes:" << node_num << endl;
    cout << "Number of Vertex:" << vertex_num << endl;
//    total_size_trie = 0;
//    total_size_trie = node_num * (4 + 4 + 8 + 8 + 1) + vertex_num * 4;
    for (int i = 0; i < trie_vertex_node.size(); i++) {
        for (int temp_k = 1; temp_k <= k_max; temp_k++) {
            total_size_trie += trie_vertex_node[i][temp_k].size() * 8;
        }
    }
    cout << "Trie Size:" << (double)total_size_trie / 1024.0 / 1024.0 << "MB" << endl;
    index_size = (double)total_size_trie / 1024.0 / 1024.0;

    long long total_graph_size = 0;
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        for (auto l = 0; l < number_of_layers; l++) {
            total_graph_size += multilayer_graph.adjacency_list[v][l].size() * 4;
        }
    }
    cout << "Graph Size:" << (double)total_graph_size / 1024.0 / 1024.0 << "MB" << endl;
    graph_size = (double)total_graph_size / 1024.0 / 1024.0;
}

void trie_query_test() {
    k = 1;
    trie_node *root_node;
    if (k <= tries.size() - 1) {
        root_node = tries[k];
    } else {
        return;
    }
    vector<vector<int>> layer_composition = {{0},
                                             {1},
                                             {2},
                                             {3},
                                             {1, 2},
                                             {0, 3},
                                             {1, 2, 3},
                                             {0, 2},
                                             {0, 2, 3},
                                             {1,8},
                                             {2,9},
                                             {5,6}};
    for (auto match_layers : layer_composition) {
        stack<pair<trie_node *, int>> stack_node;
        stack_node.emplace(root_node, 0);
        set<int> test_result;
        while (!stack_node.empty()) {
            auto top_element = stack_node.top();
            stack_node.pop();
            trie_node *node = top_element.first;
            int match_index = top_element.second;
            if (node->layer_id == match_layers[match_index]) {
                if (match_index < match_layers.size() - 1) {
                    match_index++;
                } else if (match_index == match_layers.size() - 1) {
                    // add vertices of the subtree into result;
                    queue<trie_node *> Q;
                    Q.push(node);
                    while (!Q.empty()) {
                        auto front = Q.front();
                        Q.pop();
                        for (auto v: front->vertex_vec) {
                            test_result.insert(v);
                        }
                        for (auto child_node: front->child_node) {
                            Q.push(child_node);
                        }
                    }
                    continue;
                }
            }
            for (int i = node->child_node.size() - 1; i >= 0; i--) {
                if (node->child_node[i]->layer_id <= match_layers[match_index]) {
                    stack_node.emplace(node->child_node[i], match_index);
                }
            }
        }
        cout << "trie count:" << test_result.size() << endl;

        int count_num = 0;
        for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
            if (k > PMC_index[v].size() - 1 || PMC_index[v][k].size() == 0) {
                continue;
            }
            bool exist = false;
            for (int layer_num = match_layers.size(); layer_num < number_of_layers + 1; layer_num++) {
                if (exist || layer_num > PMC_index[v][k].size() - 1)
                    break;
                for (auto layers: PMC_index[v][k][layer_num]) {
                    bool dominate = true;
                    for (auto layer: match_layers) {
                        if (std::find(layers.begin(), layers.end(), layer) == layers.end()) {
                            dominate = false;
                            break;
                        }
                    }
                    if (dominate) {
                        exist = true;
                        break;
                    }
                }
            }
            if (exist) {
                count_num++;
            }
        }
        cout << "index count:" << count_num << endl;
        if (test_result.size() != count_num)
            cout << "!!!";
    }
}

void MLCS_trie(int query_vertex) {
    if (k > k_max)
        return;
    // get the layer compositions of the query vertex
    vector<vector<int>> layer_composition;
    for (auto leaf_node : trie_vertex_node[query_vertex][k]) {
        vector<int> temp_layers;
        auto node = leaf_node;
        while (node != nullptr) {
            if (node->layer_id != -1) {
                temp_layers.push_back(node->layer_id);
            } else {
                break;
            }
            node = node->father_node;
        }
        std::reverse(temp_layers.begin(), temp_layers.end());
        if (temp_layers.size() >= d) {
            layer_composition.push_back(temp_layers);
        }
    }
    set<vector<int>> candidate_layers_set;
    for (auto layers: layer_composition) {
        auto sub_layers = generate_sub_layers(layers);
        for (auto sub_layer: sub_layers) {
            candidate_layers_set.insert(sub_layer);
        }
    }
//    set<vector<int>> candidate_layers_set;
//    for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
//        for (auto layers : PMC_index[query_vertex][k][layer_num]) {
//            auto sub_layers = generate_sub_layers(layers);
//            for (auto sub_layer : sub_layers) {
//                candidate_layers_set.insert(sub_layer);
//            }
//        }
//    }
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1);
    vector<bool> vertex_visited(multilayer_graph.maximum_vertex + 1);
    for (auto match_layers : candidate_layers_set) {
        std::fill(vertex_exist.begin(), vertex_exist.end(), false);
        std::fill(vertex_visited.begin(), vertex_visited.end(), false);

        stack<pair<trie_node *, int>> stack_node;
        trie_node* root_node = tries[k];
        stack_node.emplace(root_node, 0);
        // get all vertices in match layers
        while (!stack_node.empty()) {
            auto top_element = stack_node.top();
            stack_node.pop();
            trie_node *node = top_element.first;
            int match_index = top_element.second;
            if (node->layer_id == match_layers[match_index]) {
                if (match_index < match_layers.size() - 1) {
                    match_index++;
                } else if (match_index == match_layers.size() - 1) {
                    // add vertices of the subtree into result;
                    queue<trie_node *> Q;
                    Q.push(node);
                    while (!Q.empty()) {
                        auto front = Q.front();
                        Q.pop();
                        for (auto v: front->vertex_vec) {
                            vertex_exist[v] = true;
                        }
                        for (auto child_node: front->child_node) {
                            Q.push(child_node);
                        }
                    }
                    continue;
                }
            }
            for (int i = node->child_node.size() - 1; i >= 0; i--) {
                if (node->child_node[i]->layer_id <= match_layers[match_index]) {
                    stack_node.emplace(node->child_node[i], match_index);
                }
            }
        }

        // using BFS search from the query vertex
        vector<int> temp_result;
        queue<int> Q;
        Q.push(query_vertex);
        vertex_visited[query_vertex] = true;
        while (!Q.empty()) {
            int front = Q.front();
            Q.pop();
            temp_result.push_back(front);
            for (auto neighbor: projected_graph.adjacency_list[front]) {
                if (vertex_visited[neighbor]) {
                    continue;
                }
                vertex_visited[neighbor] = true;
                if (vertex_exist[neighbor]) {
                    Q.push(neighbor);
                }
            }
        }
        if (temp_result.size() > result.size())
            result = temp_result;
    }
}

void MLCS_trie_node_height(int query_vertex) {
    if (k > k_max)
        return;
    // get the layer compositions of the query vertex
    vector<vector<int>> layer_composition;
    for (auto leaf_node : trie_vertex_node[query_vertex][k]) {
        vector<int> temp_layers;
        auto node = leaf_node;
        while (node != nullptr) {
            if (node->layer_id != -1) {
                temp_layers.push_back(node->layer_id);
            } else {
                break;
            }
            node = node->father_node;
        }
        std::reverse(temp_layers.begin(), temp_layers.end());
        if (temp_layers.size() >= d) {
            layer_composition.push_back(temp_layers);
        }
    }
    set<vector<int>> candidate_layers_set;
    for (auto layers: layer_composition) {
        auto sub_layers = generate_sub_layers(layers);
        for (auto sub_layer: sub_layers) {
            candidate_layers_set.insert(sub_layer);
        }
    }
//    set<vector<int>> candidate_layers_set;
//    for (int layer_num = number_of_layers; layer_num >= d; layer_num--) {
//        for (auto layers : PMC_index[query_vertex][k][layer_num]) {
//            auto sub_layers = generate_sub_layers(layers);
//            for (auto sub_layer : sub_layers) {
//                candidate_layers_set.insert(sub_layer);
//            }
//        }
//    }
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1);
    vector<bool> vertex_visited(multilayer_graph.maximum_vertex + 1);
    for (auto match_layers : candidate_layers_set) {
        std::fill(vertex_exist.begin(), vertex_exist.end(), false);
        std::fill(vertex_visited.begin(), vertex_visited.end(), false);

        stack<pair<trie_node *, int>> stack_node;
        trie_node* root_node = tries[k];
        stack_node.emplace(root_node, 0);
        // get all vertices in match layers
        while (!stack_node.empty()) {
            auto top_element = stack_node.top();
            stack_node.pop();
            trie_node *node = top_element.first;
            int match_index = top_element.second;
            if (node->layer_id == match_layers[match_index]) {
                if (match_index < match_layers.size() - 1) {
                    match_index++;
                } else if (match_index == match_layers.size() - 1) {
                    // add vertices of the subtree into result;
                    queue<trie_node *> Q;
                    Q.push(node);
                    while (!Q.empty()) {
                        auto front = Q.front();
                        Q.pop();
                        for (auto v: front->vertex_vec) {
                            vertex_exist[v] = true;
                        }
                        for (auto child_node: front->child_node) {
                            Q.push(child_node);
                        }
                    }
                    continue;
                }
            }
            for (int i = node->child_node.size() - 1; i >= 0; i--) {
                if (node->child_node[i]->layer_id <= match_layers[match_index] && node->child_node[i]->height >= (match_layers.size() - match_index)) {
                    stack_node.emplace(node->child_node[i], match_index);
                }
            }
        }

        // using BFS search from the query vertex
        vector<int> temp_result;
        queue<int> Q;
        Q.push(query_vertex);
        vertex_visited[query_vertex] = true;
        while (!Q.empty()) {
            int front = Q.front();
            Q.pop();
            temp_result.push_back(front);
            for (auto neighbor: projected_graph.adjacency_list[front]) {
                if (vertex_visited[neighbor]) {
                    continue;
                }
                vertex_visited[neighbor] = true;
                if (vertex_exist[neighbor]) {
                    Q.push(neighbor);
                }
            }
        }
        if (temp_result.size() > result.size())
            result = temp_result;
    }
}

void exp() {
    ifstream test_case_file("datasets/" + dataset_name + "/" + dataset_name + "_test_case.txt");
    ofstream exp_file("datasets/" + dataset_name + "/" + dataset_name + "_index_exp.txt");
    string line;
    exp_file << "Graph Size:" << graph_size << " MB" << endl;
    exp_file << "Index Size:" << index_size << " MB" << endl;
    exp_file << "Build Time:" << build_time << " s" << endl;
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
        cout << "Start searching using index...." << endl;
        cout << "Query vertex :" << query_vertex << endl;
        exp_file << "Query vertex :" << query_vertex << endl;
        cout << "Index on " << dataset_name << endl;
        cout << "k:" << k << " " << "d:" << d << endl;
        auto start = std::chrono::high_resolution_clock::now();
//        MLCS_bitmap(query_vertex);
//        MLCS_trie(query_vertex);
        MLCS_trie_node_height(query_vertex);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "Result Num:" << result.size() << endl;
        exp_file << "Result Num:" << result.size() << endl;
        cout << "Search Time (s): " << duration.count() / 1000000000.0 <<  std::endl;
        exp_file << "Search Time (s): " << duration.count() / 1000000000.0 <<  endl << endl;
    }
    exp_file.close();
}

void get_test_case() {
    string test_case_path = "datasets/" + dataset_name + "/" + dataset_name + "_test_case.txt";
    ofstream test_case_file(test_case_path);
    for (int temp_k = 1; temp_k <= k_max; temp_k++) {
        if (tries.size() <= temp_k) {
            break;
        }
        for (int temp_d = 1; temp_d <= number_of_layers; temp_d++) {
            int query_vertex = -1;
            for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                if (PMC_index[v].size() >= temp_k +1) {
                    for (auto temp_d_d = temp_d; temp_d_d <= number_of_layers; temp_d_d++) {
                        if (PMC_index[v][temp_k][temp_d_d].size() != 0) {
                            query_vertex = v;
                            break;
                        }
                    }
                }
            }
            if (query_vertex != -1) {
                cout << temp_k << " " << temp_d << " " << query_vertex << endl;
                test_case_file << temp_k << " " << temp_d << " " << query_vertex << endl;
            }
        }
    }
    test_case_file.close();
}

void compress_trie() {
    int node_cnt = 0;
    for (int temp_k = 1; temp_k < tries.size(); temp_k++) {
        trie_node* root = tries[temp_k];
        stack<trie_node*> node_stack;
        for (auto node : root->child_node) {
            node_stack.emplace(node);
        }
        node_cnt++;
        while (!node_stack.empty()) {
            trie_node* node = node_stack.top();
            node_cnt++;
            node_stack.pop();
            bool compress = false;
            if (node->child_node.size() == 1) {
                trie_node* temp_node = node->child_node[0];
//                cout << temp_node->vertex_vec.size() << endl;
                while (temp_node->vertex_vec.empty()) {
                    node->layer_id_vec.push_back(temp_node->layer_id);
                    if (temp_node->child_node.size() == 1) {
                        temp_node = temp_node->child_node[0];
                        compress = true;
                    } else {
                        break;
                    }
                }
                if (compress) {
                    node->child_node[0] = temp_node;
                    node->compressed = true;
                    node->layer_id_vec.insert(node->layer_id_vec.begin(), node->layer_id);
                }
            }
            for (auto child_node : node->child_node) {
                node_stack.emplace(child_node);
            }
        }
    }
    cout << "Node after compressing : " << node_cnt << endl;
}

int main(int argc, char* argv[]) {
    srand(time(0));
    clock_t start_time, end_time;
    dataset_name = argv[1];
    string dataset_path = "datasets/" + dataset_name + "/" + dataset_name + ".txt";
    cout << "Loading graph...." << endl;
    MultilayerGraph multilayer_graph_dataset(dataset_path);
    multilayer_graph = multilayer_graph_dataset;
    ProjectedGraph projected_graph_dataset(multilayer_graph);
    projected_graph = projected_graph_dataset;
    number_of_layers = multilayer_graph.number_of_layers;

    cout << "Loading index...." << endl;
    load_index();
    auto start = std::chrono::high_resolution_clock::now();
    build_trie();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    build_time = duration.count() / 1000000000.0;

//    compress_trie();
//    trie_query_test();
//    cout << "Convert index to bitmap...." << endl;
//    convert_bitmap_index();
//    index_correct_test();
//    get_test_case();
    exp();

//    result_test();

//    vector<int> k_list = {1, 2, 3, 4};
//    vector<int> d_list = {1, 2, 3, 4, 5};
//    vector<int> k_list = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//    vector<int> d_list = {1, 2, 3, 4};
//    vector<int> k_list = {1};
//    vector<int> d_list = {1, 2, 3};
//    vector<pair<int, int>> k_d_list = {{5,1},{5,2},{5,3},{5,4},{5,5},{5,6},{10,1},{10,2},{10,3},{10,4},{10,5},{15,1},{15,2},{15,3},{15,4},{20,1},{20,2},{20,3}};
//
//    ofstream exp_file("datasets/" + dataset_name + "/" + dataset_name + "_index_exp.txt");
//    for (auto k_d_pair : k_d_list) {
//            result.clear();
//            k = k_d_pair.first;
//            d = k_d_pair.second;
//            exp_file << "k:" << k << " " << "d:" << d << endl;
//            cout << "Running baseline...." << endl;
//            auto start = std::chrono::high_resolution_clock::now();
//            vector<int> test_result = baseline();
//            auto end = std::chrono::high_resolution_clock::now();
//            auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//            cout << "Search Time : " << duration1.count() / 1000.0 << "s" << std::endl;
//            exp_file << "Search Time : " << duration1.count() / 1000.0 << "s" << endl << endl;
//            if (test_result.size() == 0) {
//                cout << "No result!" << endl;
//                exp_file << "No result!" << endl << endl;
//                continue;
//            }
//            cout << "Start searching using index...." << endl;
//            int query_vertex = test_result[0];
//            cout << "Query vertex :" << query_vertex << endl;
//            exp_file << "Query vertex :" << query_vertex << endl;
//            start = std::chrono::high_resolution_clock::now();
//            MLCS_bitmap(query_vertex);
//            end = std::chrono::high_resolution_clock::now();
//            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
//            cout << "Result Num:" << result.size() << endl;
//            exp_file << "Result Num:" << result.size() << endl;
//            cout << "Search Time : " << duration.count() / 1000000.0 << "ms" << std::endl;
//            exp_file << "Search Time : " << duration.count() / 1000000.0 << "ms" << endl << endl;
//    }
//    exp_file.close();

////    k = 5;
////    d = 2;
//
////    cout << "Running baseline...." << endl;
////    vector<int> test_result = baseline();
////    if (test_result.size() == 0) {
////        cout << "No result!" << endl;
////        return 0;
////    }
//
//    cout << "Start searching using index...." << endl;
//    int query_vertex = 223466;
////    int query_vertex = test_result[test_result.size()-1];
////    cout << "Query vertex :" << test_result[1] << endl;
//    cout << "Query vertex :" << query_vertex << endl;
//    auto start = std::chrono::high_resolution_clock::now();
////    MLCS(query_vertex);
//    MLCS_bitmap(query_vertex);
//    auto end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//    std::cout << "Search Time : " << duration.count() << " microseconds" << std::endl;

    return 0;
}