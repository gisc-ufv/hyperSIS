from dataclasses import dataclass
from typing import Literal, Tuple, Union, Optional

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
    network_file : str
        Path to the network edgelist file.
    network_format : str
        Format of the input network file. Options are:
        'edgelist' (default), 'bipartite', 'xgi_json', or 'hif'.
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
    network_file: str = "example.edgelist"
    network_format: Literal["edgelist", "bipartite", "xgi", "hif"] = "edgelist"
    output_dir: Optional[str] = None

    # Algorithm
    algorithm: Literal["HB_OGA", "NB_OGA"] = "HB_OGA"
    sampler: Literal["rejection_maxheap", "btree"] = "rejection_maxheap"

    # Dynamics
    tmax: int = 100
    use_qs: bool = False
    n_samples: int = 10
    time_scale: Literal["uniform", "powerlaw"] = "uniform"
    initial_condition: Union[Tuple[Literal["fraction"], float], Tuple[Literal["number"], int]] = ("fraction", 1.0)

    # IO export
    export_states: bool = False

    # Epidemic
    par_b: float = 0.5
    par_theta: float = 0.5
