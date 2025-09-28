from .core import run_simulation
from .types import SimulationArgs

__version__ = "1.2.0"

__all__ = ["run_simulation", "SimulationArgs"]

# Print a text banner when the package is imported
banner = r"""
Please cite the following paper if you use hyperSIS in your research:

    Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira.
    Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks.
    arXiv:2509.20174, 2025. DOI:10.48550/arXiv.2509.20174. https://arxiv.org/abs/2509.20174
"""

print(banner)
