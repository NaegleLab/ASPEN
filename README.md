# ASPEN

## Summary
Topology reconstruction by enumeration using branch-and-bound, as described in the 2019 publication *`Publication information to come`*. ASPEN enumeration 

### Note about parallelization
Although in some cases enumeration can be done effectively on one core, this code is set up for parallel reconstruction. Because there is some overhead to setting up the parallelization framework, this code actually runs less efficiently on one core. We recommend using at least 4 worker processes to take full advantage of the parallelization.

## Requirements

### Dependencies
ASPEN topology enumeration code runs under python 2.7 and requires the `biopython` package.
The test suite additionally requires the `mock` package.

### Input requirements
A file containing pairwise path length distributions (histograms) between each pair of leaves to be incorporated into reconstructed topologies.

Examples of path length distributions with short enumeration run times are provided.

## Basic operation instructions

1. Make sure `aspen` is in the PYTHONPATH
2. Place pairwise path length distributions file in working directory
3. Call `python -m aspen.run <`*`distributions file`*`> [`*`number of processes`*`]`

Default number of processes to use for enumeration is 4.

## License
ASPEN is published under the GNU General Public License as published by the Free Software Foundation, version 3 of the License.  See license.txt