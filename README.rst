================================================================================
Complex Adaptive Systems Assignment 2
================================================================================

This repository holds the code for our (Ariana Villegas, Christopher Leap,
Emmanuel Ohiri) report for assignment 2.

Part 1
--------------------------------------------------------------------------------
To run the code for part 1::

        $ foo bar


Part 2
--------------------------------------------------------------------------------
To run the code for part 2, first activate the conda environment::

        $ conda env create -f=./part2/environment.yml
        $ conda activate antigenic

Then run::

        $ python -m part2

This will perform a random walk through the neutral network. From each node
in the neutral network, it will generate ``breadth`` random mutations. If a
random mutation is ``max-dist`` mutations away from the neutral network, it is
marked as a leaf of the search tree and the walk is stopped down that path.
This will generate a depth-first-search tree of the neutral network, rooted at
the original genome, with a maximum depth of ``max-depth``.

For more information about the command-line arguments, run::

        $ python -m part2 --help

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The data used for fitness calculation was pulled from https://github.com/jbloomlab/SARS2-mut-fitness/
  using::

        $ curl https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/main/results/aa_fitness/aamut_fitness_all.csv > part2/aamut_fitness_all.csv

  and then converted using::

        $ python -m part2.generate_rbd_fitness

  For more information about command-line arguments, run::

        $ python -m part2.generate_rdb_fitness --help

Part 3
--------------------------------------------------------------------------------
To run the code for part 3::

        $ foo bar
