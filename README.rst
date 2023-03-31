================================================================================
Complex Adaptive Systems Assignment 2
================================================================================

This repository holds the code for our (Ariana Villegas, Christopher Leap,
Emmanuel Ohiri) report for assignment 2.

Conda Environment
--------------------------------------------------------------------------------
Before running any of the code, create and activate our conda environment::

        $ conda env create -f=environment.yml
        $ conda activate antigenic

Part 1
--------------------------------------------------------------------------------
For A section, run::

        $ python -m part1.a

This will perform an approximate count of synonymous mutations in a ``n`` 
neutral network. The process used to count the mutations is explained in 
the report.

For more information about the command-line arguments, run::

        $ python -m part1.a --help

For B section, run::

        $ python -m part1.b

This will perform an approximate count of antigenically neutral mutations 
in a ``n`` neutral network. The process used to count the neutral mutations 
is explained in the report.

For more information about the command-line arguments, run::

        $ python -m part1.b --help

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The data used for escape table was pulled from https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps
  using::

        $ curl https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv > part1/escape_calculator_data.csv

  and then converted using::

        $ python -m part1.process_escape_table

  For more information about command-line arguments, run::

        $ python -m part1.process_escape_table --help


Part 2
--------------------------------------------------------------------------------
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
