===========================================================

Improved TCH

===========================================================


Introduction
===========================================================
An implementation of Improved Time-dependent Contraction Hierarchies (TCH). 
TCH is a fast and exact solver for pathfinding in Time-dependent road networks.
For more details, please see [1]. This implementation further improve the TCH 
query performance by (i) applying landmark heuristic; (ii) applying (R)TCPD heuristic; 
and (iii) splitting the time domain to construct a TCH for different time periods.
The implementation of Improved TCH is mainly following the description by Shen et al. [2]

It must be noted that Improved TCH is released under the terms of GNU AGPL
version 3 (see files 'LICENSE' and 'COPYING'). Also, note that the implementation of Improved
TCH is based on the KaTCH (https://github.com/GVeitBatz/KaTCH). Please also see their license
when using the source code. 

Finally, the source code of Improved TCH is released mainly for research purposes. 
If you are using our source code, please carefully cite our paper [2].




Dataset
===========================================================
Improved TCH contains a real-world dataset taken from the public repository (https://uofi.app.box.com/v/NYC-traffic-estimates)
The dataset contains the road network for New York (NY) and the historical travel time that is estimated every hour 
during the entire 2013 year. To compute the travel time function for each edge, we take the travel time data from Tuesday to Thursday 
and average them for each hour after filtering out the data by two standard deviations.

The converted NY dataset contains two files: (i) NY.tpgr (input file for building TCH) and (ii) NY.coordinate (The location of each vertex).
For other synthetic datasets used in [2], please contact us if you need them.

Requirements
===========================================================

Libraries
----------------
- OpenMP
- BOOST strong typedef and pairing heap

Language Version
----------------
- C++14



Compiling and Running
===========================================================
Improved TCH is currently using Cmake to compile. You need to modify
CMakeLists.txt based on your machine setting. After that, run
"make fast" to compile.

Currently, we provide two bash scripts to quickly reproduce the
experimental results reported in paper [2].

- bash run_all_preprocessing.sh [DATASET_NAME] <br />
e.g., run "bash run_all_preprocessing.sh dataset/NY" <br />
This bash command creates all the indexes (e.g., TCH, STCH, MTCH, landmarks and (R)TCPDs)
needed for Improved TCH.

- bash run_all_experiments.sh [DATASET_NAME]  <br />
e.g., run "bash run_all_experiments.sh dataset/NY" <br />
This bash command simply runs all the algorithms described in [2] using randomly generated queries.
It will automatically output the query performance of each algorithm into the directory: dataset/NY/results/performance.
For specific algorithm, please see ./run_experiments.cpp for detail.

Contact
===========================================================
For any question, please contact Bojie.Shen@monash.edu.



References
==========

[1] Gernot Veit Batz, Robert Geisberger, Peter Sanders, and Christian
    Vetter. "Minimum Time-Dependent Travel Times with Contraction
    Hierarchies", ACM Journal of Experimental Algorithmics,
    18(1.4):1--43, 2013.
    
[2] B. Shen, M. A. Cheema, D. Harabor, P. J. Stuckey,
    Improving Time-dependent Contraction Hierarchies,
    in: Proceedings of the 32nd International Conference on 
    Automated Planning and Scheduling, ICAPS 2022.
    
