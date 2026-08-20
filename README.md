# Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks

> [!NOTE]
> Code implemented using the [Fortran Package Manager](https://fpm.fortran-lang.org/).
>
> Main paper: *Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks*, by Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira.
>
> Reference: [Nature Communications **17**, 8665 (2026)](https://www.nature.com/articles/s41467-026-75402-0) DOI:10.1038/s41467-026-75402-0
>
> Preprint: [arxiv:2509.20174](https://arxiv.org/abs/2509.20174) [DOI:10.48550/arXiv.2509.20174](https://arxiv.org/abs/2509.20174)

[![Nature Communications **17**, 8665 (2026)](natcomm.png)](https://www.nature.com/articles/s41467-026-75402-0)

## Hyper-SIS Dynamical Model

This code simulates SIS dynamics on hypergraphs (Hyper-SIS). Each of the $N$ agents can be either susceptible ($\sigma_i = 0$) or infected ($\sigma_i = 1$). Infections occur via hyperedges, which are active if a critical mass of members is infected, while infected nodes recover spontaneously.

Key points:

- Node recovery rate: $\alpha = 1$.
- Hyperedge activation threshold: $\theta(m) = 1 + (m-1)\theta_0$, where $m$ is the hyperedge order.
- Infection rate as a function of hyperedge order: $\beta(m) = \beta[1 + b(m-1)]$.
- Pairwise infection rate: $\beta(1) = \beta$.
- Parameters `par_b` and `par_theta` correspond to $b$ and $\theta_0$.

See the main paper for full details.

## Using it as a Fortran dependency

Add this package as a dependency using the [Fortran Package Manager](https://fpm.fortran-lang.org/) (fpm):

```toml
[dependencies]
hyperSIS.git = "https://github.com/gisc-ufv/hyperSIS"
```

See the [documentation](http://pages.giscbr.org/hyperSIS/) and main program for details.

## Python package

> [!TIP]
> 💡 ***A Google Colab notebook demonstrating all installation and usage steps is available [here](https://colab.research.google.com/drive/1KZanZNdr1M6bEEfw0vdAGBpkHyVrV9XH?usp=sharing).***

The easiest way to use this project is through its Python interface.

### Via `PyPI` (`pip`)

Use `pip install hyperSIS` to install it.

### Build python package

In this case, you need to clone the repository manually.

Before installing, make sure that at least one Fortran compiler is available. By default, the package assumes **GNU Fortran** (`gfortran`) installed and available in your PATH. See [Installing GFortran](https://fortran-lang.org/learn/os_setup/install_gfortran/) for help.

Steps:

1. Clone the repository and enter it:

    ```sh
    git clone https://github.com/gisc-ufv/hyperSIS.git
    cd hyperSIS
    ```

2. Activate your preferred Python environment (e.g., `venv`, `conda`, etc.):

    ```sh
    # Example with venv
    python -m venv venv
    source venv/bin/activate

    # Example with conda
    conda create -n hyperSIS python=3.11
    conda activate hyperSIS
    ```

3. Install the Python package:

    ```sh
    pip install ./python
    ```

    - If you want to use another compiler and/or Fortran flags, set the `FC` and `FFLAGS`.

    ```sh
    # Optional: customize the Fortran compiler and flags
    export FC=gfortran # default is gfortran
    export FFLAGS="-O3 -march=native -funroll-loops" # adjust optimization flags
    pip install ./python
    ```

## Usage (Python)

> [!TIP]
> 💡 ***See [examples.ipynb](https://github.com/gisc-ufv/hyperSIS/blob/main/examples.ipynb) and [Google Colab](https://colab.research.google.com/drive/1KZanZNdr1M6bEEfw0vdAGBpkHyVrV9XH?usp=sharing) for examples.***

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
  - Enable verbose output.
  - Default: `True`

- `verbose_level: str`
  - Logging level: `'info'`, `'warning'`, `'error'`, `'debug'`.
  - Default: `warning`

- `seed: int`
  - Random seed for reproducibility.
  - Default: `42`

- `remove_files: bool`
  - Remove temporary files after execution.
  - Default: `False`

- `network: NetworkFormat`
  - Network specification as a tuple. Optional parameters are in brackets:
    - `("edgelist", path, [delimiter], [comment], [cache])`
    - `("fortran-edgelist", path, [cache])`
    - `("bipartite", path, [delimiter], [comment], [cache])`
    - `("xgi", name_or_object, [cache])`
    - `("xgi_json", path, [cache])`
    - `("hif", path, [cache])`
    - `("PL", gamma, N, [sample])`
    - `("networkx", nx.Graph, [cache])`
  - Default: `("PL", 3.0, 100, 1)`

- `output_dir: Optional[str]`
  - Directory to store simulation output. If `None`, a temporary folder is used.
  - Default: `None`

- `algorithm: str`
  - Simulation algorithm: `'HB_OGA'` or `'NB_OGA'`.
  - Default: `HB_OGA`

- `sampler: str`
  - Sampling method: `'rejection_maxheap'` or `'btree'`.
  - Default: `btree`

- `tmax: int`
  - Maximum simulation time.
  - Default: `100`

- `use_qs: bool`
  - Whether to use the quasi-stationary method.
  - Default: `False`

- `n_samples: int`
  - Number of samples per simulation.
  - Default: `10`

- `time_scale: str`
  - Temporal scale for output: `'uniform'` or `'powerlaw'`.
  - Default: `uniform`

- `initial_condition: tuple`
  - Initial state specification:
    - `('fraction', float)` → fraction of infected nodes
    - `('number', int)` → exact number of initially infected nodes
  - Default: `("fraction", 1.0)`

- `export_states: bool`
  - Whether to export the full state trajectory.
  - Default: `False`

- `build_xgi_hypergraph: bool`
  - Whether to build and return the xgi hypergraph representation of the network.
  - Default: `False`

- `par_b: float`
  - Epidemic infection rate scale $b$ in $\beta(m) = \beta[1 + b(m-1)]$.
  - Default: `0.5`

- `par_theta: float`
  - Epidemic critical mass threshold $\theta_0$ in $\theta(m) = 1 + (m-1)\theta_0$.
  - Default: `0.5`

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
    - `active_states: Optional[dict]` – detailed active states per sample and time (if requested), formatted as `{sample_id: {time: {"nodes": [...], "edges": [...]}}}`.
  - `xgi_hypergraph: Optional[xgi.core.hypergraph.Hypergraph]` – representation of the structure as an `xgi` hypergraph, if generated.

## Fortran executable `hyperSIS_sampling`

Read [hyperSIS_sampling.md](https://github.com/gisc-ufv/hyperSIS/blob/main/hyperSIS_sampling.md) for instructions.

## How to Cite

> [!IMPORTANT]
> When using this package, please cite the following paper:
>
> *Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks*, by Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira (2026)
>
> Reference: [Nature Communications **17**, 8665 (2026)](https://www.nature.com/articles/s41467-026-75402-0) DOI:10.1038/s41467-026-75402-0
>
> Preprint: [arxiv:2509.20174](https://arxiv.org/abs/2509.20174) [DOI:10.48550/arXiv.2509.20174](https://arxiv.org/abs/2509.20174)

The BibTeX entry is:

```bib
@article{Maia2026_hyperSIS,
title = {Efficient {Gillespie} algorithms for spreading phenomena in large and heterogeneous higher-order networks},
    volume = {17},
    issn = {2041-1723},
    url = {https://www.nature.com/articles/s41467-026-75402-0},
    doi = {10.1038/s41467-026-75402-0},
    number = {1},
    journal = {Nature Communications},
    author = {Maia, Hugo P. and Cota, Wesley and Moreno, Yamir and Ferreira, Silvio C.},
    month = jul,
    year = {2026},
    note = {Publisher: Nature Publishing Group},
    pages = {8665},
}
```
