# Synergetic Community Search over Large Multilayer Graphs

## Introduction
This repository contains the code used in our paper: **Synergetic Community Search over Large Multilayer Graphs**.
We propose a novel core model, called Synergetic Core (SynCore). We study the problem of Synergetic Core Community Search (SynCS) over large multilayer graphs.
## Datasets
We use 9 publicly available real-world multilayer graphs, which can be obtained from [SNAP](http://snap.stanford.edu/), [KONECT](http://konect.uni-koblenz.de/networks), and [Manlio De Domenico's Website](https://manliodedomenico.com/data.php).

An example format of the input data is shown in folder datasets/test. test.txt is the original graph, test_index.txt is the Dominant Layers-based Trie of the graph, and test_test_case.txt is the query file.

## Algorithms
The following files are the codes for our proposed algorithms.
We implemented all codes the using C++, and compiled by GCC 9.3.0.
1. **MLCS_baseline.cpp** : The basic method to handle SynCS.
2. **MLCS_prune.cpp** : The Progressive Search Algorithm (PSA) for SynCS.
3. **MLCS_decomposition.cpp** : SynCore Decomposition Algorithm, and Construction Algorithm of Dominant Layers-based Trie (DLT).
4. **MLCS_index** : DLT-based Search Algorithm (DSA).
5. **MLCS_maintenance** : Index maintenance algorithm of DLT.

## Usage
1. **Compile the project with CMakeLists.txt.**

2. **DSA.** Perform SynCS using the Progressive Search Algorithm: 

For example, after the project is compiled, just type in:
```
/SynCore/cmake-build-release/MLCS_prune test
```
3. **PSA.** Perform SynCS using the DLT-based Search Algorithm:

Perform SynCore decomposition and build the index:
```
/SynCore/cmake-build-release/MLCS_decomposition test
```
Index will be stored in test_index.txt.

Perform SynCS using the DLT-based Search Algorithm:
```
/SynCore/cmake-build-release/MLCS_index test
```

## Requirements
+ cmake
+ g++
+ OpenMP