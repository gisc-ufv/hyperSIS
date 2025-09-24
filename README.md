# Hyper-SIS

Code implemented using the [Fortran Package Manager](https://fpm.fortran-lang.org/).

Main paper: *Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks*, by Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira.

## Hyper-SIS Dynamical Model

This code simulates SIS dynamics on hypergraphs (Hyper-SIS). Each of the $N$ agents can be either susceptible ($\sigma_i = 0$) or infected ($\sigma_i = 1$). Infections occur via hyperedges, which are active if a critical mass of members is infected, while infected nodes recover spontaneously.

Key points:

- Node recovery rate: $\alpha = 1$.
- Hyperedge activation threshold: $\theta(m) = 1 + (m-1)\theta_0$, where $m$ is the hyperedge order.
- Infection rate as a function of hyperedge order: $\beta(m) = \beta[1 + b(m-1)]$.
- Pairwise infection rate: $\beta(1) = \beta$.
- Parameters `par_b` and `par_theta` correspond to $b$ and $\theta_0$.

See the main paper for full details.

## Install



## Usage

*See [examples.ipynb](https://github.com/gisc-ufv/hyperSIS/blob/main/examples.ipynb) for examples.*

Import the package with

```python
import hyperSIS as hs
```

The simulation interface revolves around **two main objects**:

1. `SimulationArgs`
   A dataclass containing all parameters required to configure a hyperSIS simulation, including network specification, algorithm choices, temporal settings, initial conditions, and epidemic parameters.

2. `run_simulation(beta1: float, args: SimulationArgs)`
   The function that executes the simulation with the given arguments. Returns a `SimulationResult` object containing the processed results, including network mapping, temporal evolution, and statistics of infected nodes.

### Simulation arguments

The `SimulationArgs` dataclass contains all configurable parameters for running a hyperSIS simulation.

- `verbose: bool`
  Enable verbose output.
  Default: `True`

- `verbose_level: str`
  Logging level: `'info'`, `'warning'`, `'error'`, `'debug'`.
  Default: `warning`

- `seed: int`
  Random seed for reproducibility.
  Default: `42`

- `remove_files: bool`
  Remove temporary files after execution.
  Default: `False`

- `network: NetworkFormat`
  Network specification as a tuple. Optional parameters are in brackets:
  - `("edgelist", path, [delimiter], [comment], [cache])`
  - `("fortran-edgelist", path, [cache])`
  - `("bipartite", path, [delimiter], [comment], [cache])`
  - `("xgi", name_or_object, [cache])`
  - `("xgi_json", path, [cache])`
  - `("hif", path, [cache])`
  - `("PL", gamma, N, [sample])`
  Default: `("edgelist", "example.edgelist", None, "#", False)`

- `output_dir: Optional[str]`
  Directory to store simulation output. If `None`, a temporary folder is used.
  Default: `None`

- `algorithm: str`
  Simulation algorithm: `'HB_OGA'` or `'NB_OGA'`.
  Default: `HB_OGA`

- `sampler: str`
  Sampling method: `'rejection_maxheap'` or `'btree'`.
  Default: `btree`

- `tmax: int`
  Maximum simulation time.
  Default: `100`

- `use_qs: bool`
  Whether to use the quasi-stationary method.
  Default: `False`

- `n_samples: int`
  Number of samples per simulation.
  Default: `10`

- `time_scale: str`
  Temporal scale for output: `'uniform'` or `'powerlaw'`.
  Default: `uniform`

- `initial_condition: tuple`
  Initial state specification:
  - `('fraction', float)` → fraction of infected nodes
  - `('number', int)` → exact number of initially infected nodes
  Default: `("fraction", 1.0)`

- `export_states: bool`
  Whether to export the full state trajectory.
  Default: `False`

- `par_b: float`
  Epidemic infection rate scale $b$ in $\beta(m) = \beta[1 + b(m-1)]$.
  Default: `0.5`

- `par_theta: float`
  Epidemic critical mass threshold $\theta_0$ in $\theta(m) = 1 + (m-1)\theta_0$.
  Default: `0.5`

### Function

```python
run_simulation(beta1: float, args: SimulationArgs)
```

Runs a Hyper-SIS simulation on the specified network.

**Parameters:**

- `beta1: float`
  Base infection rate $\beta(1)$ for pairwise interactions.
- `args: SimulationArgs`
  Simulation parameters, including network specification, algorithm choice, number of samples, initial condition, and epidemic parameters `par_b` and `par_theta`.

**Returns:**

- `SimulationResult`
  Object containing:

  - `network: NetworkFormat` – the network specification used.
  - `node_map: dict` – mapping from original node IDs to Fortran node IDs.
  - `temporal: TemporalResult` – temporal dynamics with:
    - `t: np.ndarray` – mean time per Gillespie tick.
    - `rho_avg: np.ndarray` – mean number of infected nodes over all runs.
    - `rho_var: np.ndarray` – variance of infected nodes.
    - `n_samples: int` – number of runs where infection is non-zero.

## How to Cite

When using this package, please cite the following paper:

*Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks*, by Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira (2025)

The BibTeX entry is:

[to be defined]
