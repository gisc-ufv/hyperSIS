from dataclasses import dataclass
from typing import Literal, Tuple, Union, Optional, Dict, Any
import numpy as np
import xgi

NetworkFileResult = Tuple[str, Dict[Any, int]]

NetworkFormatEdgelist = Tuple[
    Literal["edgelist"], str, Optional[str], Optional[str], Optional[bool]
]
# ("edgelist", path, delimiter, comment, cache)

NetworkFormatFortranEdgelist = Tuple[
    Literal["fortran-edgelist"], str
]
# ("fortran-edgelist", path, cache)

NetworkFormatBipartite = Tuple[
    Literal["bipartite"], str, Optional[str], Optional[str], Optional[bool]
]
# ("bipartite", path, delimiter, comment, cache)

NetworkFormatXGI = Tuple[
    Literal["xgi"], Union[str, xgi.core.hypergraph.Hypergraph]
]
# ("xgi", name_or_object)

NetworkFormatJSON = Tuple[
    Literal["xgi_json"], str, Optional[bool]
]
# ("xgi_json", path, cache)

NetworkFormatHIF = Tuple[
    Literal["hif"], str, Optional[bool]
]
# ("hif", path, cache)

NetworkFormatPL = Tuple[
    Literal["PL"],
    Literal[2.5, 3.0, 6.0],
    Literal[100, 1000, 10000, 100000, 1000000],
    Optional[Literal[1, 2, 3, 4, 5]]
]
# ("PL", gamma, N, sample)

NetworkFormat = Union[
    NetworkFormatEdgelist,
    NetworkFormatFortranEdgelist,
    NetworkFormatBipartite,
    NetworkFormatXGI,
    NetworkFormatJSON,
    NetworkFormatHIF,
    NetworkFormatPL
]

@dataclass
class TemporalResult:
    t: np.ndarray
    rho_avg: np.ndarray
    rho_var: np.ndarray
    n_samples: int

@dataclass
class SimulationResult:
    network: NetworkFormat
    node_map: dict
    temporal: TemporalResult
    active_states: Optional[dict] = None  # {sample_id: {time: {"nodes": [...], "edges": [...]}}}
    xgi_hypergraph: Optional[xgi.core.hypergraph.Hypergraph] = None

@dataclass
class SimulationArgs:
    """
    Simulation arguments for hyperSIS

    Attributes
    ----------
    verbose : bool
        Enable verbose output.
    verbose_level : str
        Logging level: 'info', 'warning', 'error', 'debug'.
    seed : int
        Random seed.
    remove_files : bool
        Remove temporary files after execution.
    network: NetworkFormat
        Network specification as a tuple ([optional parameters in brackets]):
        ("edgelist", path, [delimiter], [comment], [cache])
        ("fortran-edgelist", path, [cache])
        ("bipartite", path, [delimiter], [comment], [cache])
        ("xgi", name_or_object, [cache])
        ("xgi_json", path, [cache])
        ("hif", path, [cache])
        ("PL", gamma, N, [sample])
    output_dir : Optional[str]
        Directory to store simulation output. If None, uses temporary folder.
    algorithm : str
        Simulation algorithm: 'HB_OGA' or 'NB_OGA'.
    sampler : str
        Sampling method: 'rejection_maxheap' or 'btree'.
    tmax : int
        Maximum simulation time.
    use_qs : bool
        Whether to use quasi-stationary method.
    n_samples : int
        Number of samples per simulation.
    time_scale : str
        Temporal scale for output: 'uniform' or 'powerlaw'.
    initial_condition : tuple
        Initial state specification:
        ('fraction', float) -> fraction of infected nodes
        ('number', int) -> exact number of initially infected nodes
    export_states : bool
        Whether to export the full state trajectory.
    par_b : float
        Epidemic parameter beta.
    par_theta : float
        Epidemic parameter theta.
    """
    # General
    verbose: bool = True
    verbose_level: Literal["info", "warning", "error", "debug"] = "warning"
    seed: int = 42
    remove_files: bool = False

    # IO
    network: NetworkFormat = ("PL", # default
                              3.0,
                              "100",
                              None
                            )
    output_dir: Optional[str] = None

    # Algorithm
    algorithm: Literal["HB_OGA", "NB_OGA"] = "HB_OGA"
    sampler: Literal["rejection_maxheap", "btree"] = "btree"

    # Dynamics
    tmax: int = 100
    use_qs: bool = False
    n_samples: int = 10
    time_scale: Literal["uniform", "powerlaw"] = "uniform"
    initial_condition: Union[Tuple[Literal["fraction"], float], Tuple[Literal["number"], int]] = ("fraction", 1.0)

    # IO export
    export_states: bool = False
    build_xgi_hypergraph: bool = False

    # Epidemic
    par_b: float = 0.5
    par_theta: float = 0.5
