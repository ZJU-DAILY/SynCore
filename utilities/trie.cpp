//
// Created by 罗程阳 on 2024/4/18.
//

#include "trie.h"

trie_node::trie_node() {

}

trie_node::trie_node(int id) {
    layer_id = id;
    father_node = nullptr;
    height = 1;
}
