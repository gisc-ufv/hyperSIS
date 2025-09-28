# Run temporal dynamics on a hypergraph network

Manual page for `Run temporal dynamics on a hypergraph network` version 1.0

`Run temporal dynamics on a hypergraph network [--output value] [--remove-files value] [--edges-file value] [--algorithm value] [--sampler value] [--tmax value] [--use-qs value] [--n-samples value] [--time-scale value] [--initial-fraction value] [--initial-number value] --beta1 value [--par-b value] [--par-theta value] [--export-states value] [--seed value] [--verbose value] [--verbose-level value] [--help] [--markdown] [--version]`

Sep 2025

### Short description

This program runs a temporal dynamics on a hypergraph network, using the Gillespie algorithm. The network can be read from a file or generated as a small example. The dynamics can be configured with various parameters.

### Command line options:


Required switches:

* `--beta1 value`, `-b1 value`
    Infection rate parameter beta1

Optional switches:

* `--output value`, `-o value`
    default value ./output
    Output prefix for result files

* `--remove-files value`, `-rm value`
    default value false
    Remove existing output files before running

* `--edges-file value`, `-e value`
    File containing the edge list of the network as input

* `--algorithm value`, `-a value`  , value in: `HB_OGA,NB_OGA`
    default value HB_OGA
    Dynamics algorithm to use

* `--sampler value`, `-s value`  , value in: `rejection_maxheap,btree`
    default value rejection_maxheap
    Sampler choice

* `--tmax value`, `-t value`
    default value 1000.0
    Maximum simulation time

* `--use-qs value`, `-qs value`
    default value false
    Use Quasi-Stationary method

* `--n-samples value`, `-ns value`
    default value 10
    Number of samples to average over

* `--time-scale value`, `-ts value`  , value in: `powerlaw,uniform`
    default value powerlaw
    Time scale to use

* `--initial-fraction value`, `-if value`
    default value 1.0
    Initial fraction of infected nodes

* `--initial-number value`, `-in value`
    default value 0
    Initial number of infected nodes (overrides initial-fraction)

* `--par-b value`, `-pb value`
    default value 0.5
    Dynamical parameter b

* `--par-theta value`, `-pt value`
    default value 0.5
    Dynamical parameter theta

* `--export-states value`, `-es value`
    default value false
    Export the states of nodes and edges at the end of each sample (be careful with large networks)

* `--seed value`, `-rs value`
    default value 42
    Random seed for the dynamics

* `--verbose value`, `-vv value`
    default value true
    Enable verbose logging

* `--verbose-level value`, `-vl value`  , value in: `error,warning,info,debug`
    default value info
    Logging level

* `--help`, `-h`
    Print this help message

* `--markdown`, `-md`
    Save this help message in a Markdown file

* `--version`, `-v`
    Print version
