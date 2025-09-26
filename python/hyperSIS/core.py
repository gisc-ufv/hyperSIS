from pathlib import Path
import subprocess
import tempfile
from typing import Tuple
import numpy as np

from .types import SimulationArgs, TemporalResult, SimulationResult
from .io_utils import process_results, prepare_network_file

def run_simulation(beta1: float, args: SimulationArgs) -> SimulationResult:
    """
    Runs a simulation by calling the Fortran binary via fpm.

    Parameters
    ----------
    beta1 : float
        Value of the $\beta(1)$ parameter.
    args : SimulationArgs
        Simulation arguments.

    Returns
    -------
    times : np.ndarray
        Array of time points.
    rho_mean : np.ndarray
        Mean infection density.
    """
    # If no output_dir is provided, create a temporary one
    if args.output_dir is None:
        tmp_context = tempfile.TemporaryDirectory()
        tmpdir = Path(tmp_context.name)
        cleanup_tmp = True
    else:
        tmpdir = Path(args.output_dir)
        tmpdir.mkdir(parents=True, exist_ok=True)
        cleanup_tmp = False

    try:
        # Here we prepare the network file for Fortran
        network_file_fortran, map_nodes = prepare_network_file(
            args.network
        )

        xgi_hypergraph = None
        if (args.build_xgi_hypergraph):
            # based on the network_file_fortran, build the xgi hypergraph and store it in the result
            # use xgi to read edgelist and assume node_type as int
            import xgi
            xgi_hypergraph = xgi.read_edgelist(network_file_fortran, nodetype=int)

        # Chooses initial_fraction or initial_number
        kind, value = args.initial_condition
        if kind == "number":
            initial_arg = f"--initial-number {value}"
        else:
            initial_arg = f"--initial-fraction {value}"

        # Builds the command string
        inner_cmd = (
            f"fpm run hyperSIS_sampling --compiler='{args.COMPILER}' --flag='{args.COMPILER_FLAGS}' -- "
            f"--output {tmpdir}/ "
            f"--remove-files {args.remove_files} "
            f"--edges-file {network_file_fortran} "
            f"--algorithm {args.algorithm} "
            f"--sampler {args.sampler} "
            f"--tmax {args.tmax} "
            f"--use-qs {args.use_qs} "
            f"--n-samples {args.n_samples} "
            f"--time-scale {args.time_scale} "
            f"{initial_arg} "
            f"--beta1 {beta1} "
            f"--par-b {args.par_b} "
            f"--par-theta {args.par_theta} "
            f"--export-states {args.export_states} "
            f"--seed {args.seed} "
            f"--verbose {args.verbose} "
            f"--verbose-level {args.verbose_level} "
        )

        # Executes the command
        result = subprocess.run(inner_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Command failed with code {result.returncode}\n"
                f"stderr:\n{result.stderr}"
            )

        # Processes the results and returns numpy arrays
        times, rho_mean, rho_var, n_samples, active_states = process_results(tmpdir, args.export_states)

        # Create TemporalResult dataclass instance
        temporal = TemporalResult(
            t=times,          # mean time per Gillespie tick
            rho_avg=rho_mean,  # mean number of infected nodes (over all runs)
            rho_var=rho_var,   # variance of infected nodes
            n_samples=n_samples  # number of runs where rho != 0
        )

        # Create SimulationResult dataclass instance
        result = SimulationResult(
            network=args.network, # original network specification
            node_map=map_nodes, # mapping from original node IDs to Fortran node IDs
            temporal=temporal,
            active_states=active_states,  # full state trajectories if export_states is True
            xgi_hypergraph=xgi_hypergraph  # XGI hypergraph representation if built
        )

        # Return the structured result object
        return result

    finally:
        if cleanup_tmp:
            tmp_context.cleanup()
