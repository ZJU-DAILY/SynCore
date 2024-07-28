//
// Created by 罗程阳 on 2024/4/28.
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
string dataset_name = "homo";
//string dataset_name = "sacchcere";
//string dataset_name = "stack";

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

using namespace std;

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

void build_trie() {
    tries.clear();
    tries.resize(k_max + 1);
    trie_vertex_node.clear();
    trie_vertex_node.resize(multilayer_graph.maximum_vertex + 1);
    for (int i = 0; i < multilayer_graph.maximum_vertex + 1; i++) {
        trie_vertex_node[i].resize(k_max + 1);
    }
    long long total_size_trie = 0;
    int node_num = 0;
    int vertex_num = 0;
    for (int temp_k = 1; temp_k <= k_max; temp_k++) {
        auto *root = new trie_node(-1);
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
//    cout << "Number of Trie Nodes:" << node_num << endl;
//    cout << "Number of Vertex:" << vertex_num << endl;
//    total_size_trie = 0;
//    total_size_trie = node_num * (4 + 8 + 8 + 1) + vertex_num * 4;
    for (int i = 0; i < trie_vertex_node.size(); i++) {
        for (int temp_k = 1; temp_k <= k_max; temp_k++) {
            total_size_trie += trie_vertex_node[i][temp_k].size() * 8;
        }
    }
//    cout << "Trie Size:" << (double)total_size_trie / 1024.0 / 1024.0 << "MB" << endl;
}

void insert_test(ofstream &exp_file) {
    vector<int> num_list = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
//    vector<string> num_list = {"0.1", "0.2", "0.3", "0.4", "0.5"};
    for (auto num : num_list) {
        vector<vector<vector<vector<vector<int>>>>> temp_PMC_index;
        temp_PMC_index.resize(multilayer_graph.maximum_vertex + 1);
        ifstream index_file("datasets/" + dataset_name + "/" + dataset_name + "_index_insert_" + to_string(num) + ".txt");
        string line;
        int current_v = -1;
        int current_k = -1;
        int new_k_max = -1;
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
                temp_PMC_index[current_v].push_back({});
            } else if ((words[0] == "k")) {
                current_k = stoi(words[1]);
                new_k_max = max(k_max, current_k);
                temp_PMC_index[current_v].push_back({});
            } else {
                vector<int> layers;
                for (auto str : words) {
                    layers.push_back(stoi(str));
                }
                if (temp_PMC_index[current_v][current_k].size() == 0) {
                    temp_PMC_index[current_v][current_k].resize(number_of_layers + 1);
                }
                temp_PMC_index[current_v][current_k][layers.size()].push_back(layers);
            }
        }
        index_file.close();

        vector<map<int, vector<vector<int>>>> k_updated_vertex_add(new_k_max + 1);
        vector<map<int, vector<vector<int>>>> k_updated_vertex_delete(new_k_max + 1);
        vector<int> updated_vertex;
        // get the updated dominate layers
        int change_cnt = 0;
        for (int temp_k = 1; temp_k <= new_k_max; temp_k++) {
            for (auto v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                set<vector<int>> old_layers_set;
                set<vector<int>> new_layers_set;
                if (temp_k < PMC_index[v].size()) {
                    for (int l = 1; l < PMC_index[v][temp_k].size(); l++) {
                        for (const auto &layers: PMC_index[v][temp_k][l]) {
                            old_layers_set.insert(layers);
                        }
                    }
                }
                if (temp_k < temp_PMC_index[v].size()) {
                    for (int l = 1; l < temp_PMC_index[v][temp_k].size(); l++) {
                        for (const auto &layers: temp_PMC_index[v][temp_k][l]) {
                            new_layers_set.insert(layers);
                        }
                    }
                }
                vector<vector<int>> delete_layers;
                vector<vector<int>> add_layers;
                set_difference(old_layers_set.begin(), old_layers_set.end(), new_layers_set.begin(),
                               new_layers_set.end(),
                               inserter(delete_layers, delete_layers.begin()));
                set_difference(new_layers_set.begin(), new_layers_set.end(), old_layers_set.begin(),
                               old_layers_set.end(),
                               inserter(add_layers, add_layers.begin()));
                change_cnt += delete_layers.size();
                change_cnt += add_layers.size();
                if (!delete_layers.empty() || !add_layers.empty()) {
                    k_updated_vertex_delete[temp_k][v] = delete_layers;
                    k_updated_vertex_add[temp_k][v] = add_layers;
                    updated_vertex.push_back(v);
                }
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        // update the trie
        for (int temp_k = 1; temp_k <= new_k_max; temp_k++) {
            for (auto v: updated_vertex) {
                for (auto add_layers: k_updated_vertex_add[temp_k][v]) {
                    trie_node *node = tries[temp_k];
                    for (int i = 0; i < add_layers.size(); i++) {
                        int layer = add_layers[i];
                        bool exist = false;
                        for (auto child: node->child_node) {
                            if (child->layer_id == layer) {
                                node = child;
                                exist = true;
                                break;
                            }
                        }
                        if (!exist) {
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
                            new_node->father_node = node;
                            node = new_node;
                            // update the height
                            trie_node *temp_node = node;
                            while (temp_node != tries[temp_k]) {
                                if (temp_node->father_node->height < temp_node->height + 1) {
                                    temp_node->father_node->height = temp_node->height + 1;
                                    temp_node = temp_node->father_node;
                                } else {
                                    break;
                                }
                            }
                        }
                        if (i == add_layers.size() - 1) {
                            trie_vertex_node[v][temp_k].push_back(node);
                            node->vertex_vec.push_back(v);
                            node->dominate = true;
                        }
                    }
                }

                for (auto delete_layers : k_updated_vertex_delete[temp_k][v]) {
                    trie_node *node = tries[temp_k];
                    for (int layer : delete_layers) {
                        for (auto child: node->child_node) {
                            if (child->layer_id == layer) {
                                node = child;
                                break;
                            }
                        }
                    }
                    std::remove(node->vertex_vec.begin(), node->vertex_vec.end(), v);
                    std::remove(trie_vertex_node[v][temp_k].begin(), trie_vertex_node[v][temp_k].end(), node);
                    if (node->vertex_vec.empty() && node->child_node.empty()) {
                        std::remove(node->father_node->child_node.begin(), node->father_node->child_node.end(), node);
                        free(node);
                        auto temp_node = node->father_node;
                        while (temp_node != tries[temp_k]) {
                            int new_height = 0;
                            for (auto child_node : temp_node->child_node) {
                                new_height = max(new_height, child_node->height + 1);
                            }
                            if (new_height < temp_node->height) {
                                temp_node->height = new_height;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "Insert Edges: " <<  num << endl;
        cout << "Change Dominated Layers: " << change_cnt << endl;
        cout << "Updating Time (s): " << duration.count() / 1000000000.0 <<  std::endl;

        exp_file << "Insert Edges: " <<  num << endl;
        exp_file << "Change Dominated Layers: " << change_cnt << endl;
        exp_file << "Updating Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl;

        // rebuild the updated trie
        start = std::chrono::high_resolution_clock::now();
        tries.clear();
        tries.resize(k_max + 1);
        trie_vertex_node.clear();
        trie_vertex_node.resize(multilayer_graph.maximum_vertex + 1);
        for (int i = 0; i < multilayer_graph.maximum_vertex + 1; i++) {
            trie_vertex_node[i].resize(k_max + 1);
        }
        for (int temp_k = 1; temp_k <= new_k_max; temp_k++) {
            auto *root = new trie_node(-1);
            for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                if (temp_PMC_index[v].size() < temp_k + 1)
                    continue;
                for (int temp_d = 1; temp_d < temp_PMC_index[v][temp_k].size(); temp_d++) {
                    for (auto layers : temp_PMC_index[v][temp_k][temp_d]) {
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
                                new_node->father_node = node;
                                node = new_node;
                            }
                            if (i == layers.size() - 1) {
                                trie_vertex_node[v][temp_k].push_back(node);
                                node->vertex_vec.push_back(v);
                                node->dominate = true;
                            }
                        }
                    }
                }
            }
            tries[temp_k] = root;
        }
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "Rebuild Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl << endl;
        exp_file << "Rebuild Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl << endl;

        // rebuild the origin trie
        build_trie();
    }
}

void delete_test(ofstream &exp_file) {
    vector<int> num_list = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
//    vector<string> num_list = {"0.1", "0.2", "0.3", "0.4", "0.5"};
    for (auto num : num_list) {
        vector<vector<vector<vector<vector<int>>>>> temp_PMC_index;
        temp_PMC_index.resize(multilayer_graph.maximum_vertex + 1);
        ifstream index_file("datasets/" + dataset_name + "/" + dataset_name + "_index_delete_" + to_string(num) + ".txt");
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
                temp_PMC_index[current_v].push_back({});
            } else if ((words[0] == "k")) {
                current_k = stoi(words[1]);
                temp_PMC_index[current_v].push_back({});
            } else {
                vector<int> layers;
                for (auto str : words) {
                    layers.push_back(stoi(str));
                }
                if (temp_PMC_index[current_v][current_k].size() == 0) {
                    temp_PMC_index[current_v][current_k].resize(number_of_layers + 1);
                }
                temp_PMC_index[current_v][current_k][layers.size()].push_back(layers);
            }
        }
        index_file.close();

        vector<map<int, vector<vector<int>>>> k_updated_vertex_add(k_max + 1);
        vector<map<int, vector<vector<int>>>> k_updated_vertex_delete(k_max + 1);
        vector<int> updated_vertex;
        // get the updated dominate layers
        int change_cnt = 0;
        for (int temp_k = 1; temp_k <= k_max; temp_k++) {
            for (auto v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                set<vector<int>> old_layers_set;
                set<vector<int>> new_layers_set;
                if (temp_k < PMC_index[v].size()) {
                    for (int l = 1; l < PMC_index[v][temp_k].size(); l++) {
                        for (const auto &layers: PMC_index[v][temp_k][l]) {
                            old_layers_set.insert(layers);
                        }
                    }
                }
                if (temp_k < temp_PMC_index[v].size()) {
                    for (int l = 1; l < temp_PMC_index[v][temp_k].size(); l++) {
                        for (const auto &layers: temp_PMC_index[v][temp_k][l]) {
                            new_layers_set.insert(layers);
                        }
                    }
                }
                vector<vector<int>> delete_layers;
                vector<vector<int>> add_layers;
                set_difference(old_layers_set.begin(), old_layers_set.end(), new_layers_set.begin(),
                               new_layers_set.end(),
                               inserter(delete_layers, delete_layers.begin()));
                set_difference(new_layers_set.begin(), new_layers_set.end(), old_layers_set.begin(),
                               old_layers_set.end(),
                               inserter(add_layers, add_layers.begin()));
                change_cnt += delete_layers.size();
                change_cnt += add_layers.size();
                if (!delete_layers.empty() || !add_layers.empty()) {
                    k_updated_vertex_delete[temp_k][v] = delete_layers;
                    k_updated_vertex_add[temp_k][v] = add_layers;
                    updated_vertex.push_back(v);
                }
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        // update the trie
        for (int temp_k = 1; temp_k <= k_max; temp_k++) {
            for (auto v: updated_vertex) {
                for (auto add_layers: k_updated_vertex_add[temp_k][v]) {
                    trie_node *node = tries[temp_k];
                    for (int i = 0; i < add_layers.size(); i++) {
                        int layer = add_layers[i];
                        bool exist = false;
                        for (auto child: node->child_node) {
                            if (child->layer_id == layer) {
                                node = child;
                                exist = true;
                                break;
                            }
                        }
                        if (!exist) {
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
                            new_node->father_node = node;
                            node = new_node;
                            // update the height
                            trie_node *temp_node = node;
                            while (temp_node != tries[temp_k]) {
                                if (temp_node->father_node->height < temp_node->height + 1) {
                                    temp_node->father_node->height = temp_node->height + 1;
                                    temp_node = temp_node->father_node;
                                } else {
                                    break;
                                }
                            }
                        }
                        if (i == add_layers.size() - 1) {
                            trie_vertex_node[v][temp_k].push_back(node);
                            node->vertex_vec.push_back(v);
                            node->dominate = true;
                        }
                    }
                }

                for (auto delete_layers : k_updated_vertex_delete[temp_k][v]) {
                    trie_node *node = tries[temp_k];
                    for (int layer : delete_layers) {
                        for (auto child: node->child_node) {
                            if (child->layer_id == layer) {
                                node = child;
                                break;
                            }
                        }
                    }
                    std::remove(node->vertex_vec.begin(), node->vertex_vec.end(), v);
                    std::remove(trie_vertex_node[v][temp_k].begin(), trie_vertex_node[v][temp_k].end(), node);
                    if (node->vertex_vec.empty() && node->child_node.empty()) {
                        std::remove(node->father_node->child_node.begin(), node->father_node->child_node.end(), node);
                        free(node);
                        auto temp_node = node->father_node;
                        while (temp_node != tries[temp_k]) {
                            int new_height = 0;
                            for (auto child_node : temp_node->child_node) {
                                new_height = max(new_height, child_node->height + 1);
                            }
                            if (new_height < temp_node->height) {
                                temp_node->height = new_height;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "Delete Edges: " <<  num << endl;
        cout << "Change Dominated Layers: " << change_cnt << endl;
        cout << "Updating Time (s): " << duration.count() / 1000000000.0 <<  std::endl;
        exp_file << "Delete Edges: " <<  num << endl;
        exp_file << "Change Dominated Layers: " << change_cnt << endl;
        exp_file << "Updating Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl;

        // rebuild the updated trie
        start = std::chrono::high_resolution_clock::now();
        tries.clear();
        tries.resize(k_max + 1);
        trie_vertex_node.clear();
        trie_vertex_node.resize(multilayer_graph.maximum_vertex + 1);
        for (int i = 0; i < multilayer_graph.maximum_vertex + 1; i++) {
            trie_vertex_node[i].resize(k_max + 1);
        }
        for (int temp_k = 1; temp_k <= k_max; temp_k++) {
            auto *root = new trie_node(-1);
            for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
                if (temp_PMC_index[v].size() < temp_k + 1)
                    continue;
                for (int temp_d = 1; temp_d < temp_PMC_index[v][temp_k].size(); temp_d++) {
                    for (auto layers : temp_PMC_index[v][temp_k][temp_d]) {
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
                                new_node->father_node = node;
                                node = new_node;
                            }
                            if (i == layers.size() - 1) {
                                trie_vertex_node[v][temp_k].push_back(node);
                                node->vertex_vec.push_back(v);
                                node->dominate = true;
                            }
                        }
                    }
                }
            }
            tries[temp_k] = root;
        }
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        cout << "Rebuild Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl << endl;
        exp_file << "Rebuild Time (s): " << endl << duration.count() / 1000000000.0 <<  std::endl << endl;

        // rebuild the trie
        build_trie();
    }
}

int main() {
    srand(time(0));
    clock_t start_time, end_time;
    string dataset_path = "datasets/" + dataset_name + "/" + dataset_name + ".txt";
    cout << "Loading graph...." << endl;
    MultilayerGraph multilayer_graph_dataset(dataset_path);
    multilayer_graph = multilayer_graph_dataset;
    ProjectedGraph projected_graph_dataset(multilayer_graph);
    projected_graph = projected_graph_dataset;
    number_of_layers = multilayer_graph.number_of_layers;

    cout << "Loading index...." << endl;
    load_index();

    ofstream exp_file("datasets/" + dataset_name + "/" + dataset_name + "_index_maintenance.txt");
    auto start = std::chrono::high_resolution_clock::now();
    build_trie();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Building Time (s): " << duration.count() / 1000000000.0 <<  std::endl;
    exp_file << "Building Time (s): " << duration.count() / 1000000000.0 <<  std::endl;

    insert_test(exp_file);
    delete_test(exp_file);
    exp_file.close();
    return 0;
}