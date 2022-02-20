===========================================================

# Improving TCH

===========================================================


Introduction
===========================================================
An implementation of End Point Search (EPS). EPS is a fast and
exact path planner that works in Euclidean space. The implementation
of EPS is mainly follows the description by Shen et al. [1]

It must be noted that EPS is released under the terms of GNU AGPL
version 3 (see files 'LICENSE' and 'COPYING'). The source code of EPS
is released mainly for research purpose. If you are using our source
code, please carefully cite our paper [1].




Dataset
===========================================================
EPS contains four benchmark suites (dao, da2, sc1, bgmaps) retrieved
from MovingAI (https://movingai.com). The
merged-meshes are provide by the authors of Polyanya [2] and available
from repository (https://bitbucket.org/%7B3c286763-d509-45c2-b036-75814ce8955f%7D/)



Requirements
===========================================================

Libraries
----------------
- OpenMP

Language Version
----------------
- C++14



Compiling and Running
===========================================================
EPS is currently using Cmake to compile, you need to modify
CMakeLists.txt based on your machine setting. After that, run
"make fast" to compile.


Currently, we provide three bash scripts to quickly reproduce the
experimental results reported in paper [1].

- bash preprocessing.sh [MAP_NAME] <br />
e.g., run "bash preprocessing.sh dao"
This bash command creates all the indexes (visibility graph, CPD)
needed for End Point Search for all the maps in the benchmark suite (dao).

- bash benchmark_EPS.sh [MAP_NAME] <br />
e.g., run "bash benchmark_EPS.sh dao"
This bash command simply runs EPS [1] and Polyanya [2] for all the maps in
the benchmark suite (dao) using the queries taken from MovingAI.

- bash clean_index.sh <br />
e.g., run "bash clean_index.sh".
This bash command simply  delete all the indexes for all benchmark suites



Contact
===========================================================
For any question, please contact Bojie.Shen@monash.edu.



References
==========

[1] B. Shen, M. A. Cheema, D. Harabor, P. J. Stuckey,
    Euclidean pathfinding with compressed path databases,
    in: Proceedings of the Twenty-Ninth International Joint
    Conference on Artificial Intelligence, IJCAI 2020, 2020,
    pp. 4229–4235.

[2] M. Cui, D. D. Harabor, A. Grastien, Compromise-free pathfinding
    on a navigation mesh, in: Proceedings of the Twenty-Sixth International
    Joint Conference on Articial Intelligence, IJCAI 2017, 2017,
    pp. 496-502.
