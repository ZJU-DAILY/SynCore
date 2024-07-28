//
// Created by 罗程阳 on 2024/4/18.
//

#ifndef MLCS_TRIE_H
#define MLCS_TRIE_H

#include <vector>

using namespace std;
class trie_node {
public:
    int layer_id;
    vector<int> layer_id_vec;
    bool dominate = false;
    bool compressed = false;
    vector<trie_node*> child_node;
    vector<int> vertex_vec;
    trie_node* father_node;
    int height;
public:
    trie_node();
    trie_node(int id);
};
#endif //MLCS_TRIE_H
