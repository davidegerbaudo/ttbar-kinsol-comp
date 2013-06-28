ttbar-kinsol-comp
=================

Comparison of kinematic solver for ttbar

quickstart
----------

    ./runSolver.py data/bveto/105200_DIL_optimSRjets.root 
    number of solutions Betchart     :  defaultdict(<type 'int'>, {0: 15166, 2: 22888})
    time :  0.00223185352967
    number of solutions Sonnenschein :  defaultdict(<type 'int'>, {0: 27555, 2: 9543, 4: 956})
    time :  1.55088109214e-05

dileptonSolver
--------------
Compile as standalone:

    cmake .
    make

Create shared library:

    root -l -x -q makeLib.C
