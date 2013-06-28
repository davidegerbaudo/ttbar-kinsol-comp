ttbar-kinsol-comp
=================

Comparison of kinematic solver for ttbar

quickstart
----------

    ./runSolver.py data/bveto/105200_DIL_optimSRjets.root 

dileptonSolver
--------------
Compile as standalone:

    cmake .
    make

Create shared library:

    root -l -x -q makeLib.C
