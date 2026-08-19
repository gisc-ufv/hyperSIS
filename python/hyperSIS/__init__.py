from .core import run_simulation
from .types import SimulationArgs

__version__ = "1.2.3"

__all__ = ["run_simulation", "SimulationArgs"]

# Print a text banner when the package is imported
banner = r"""
Please cite the following paper if you use hyperSIS in your research:

    Hugo P. Maia, Wesley Cota, Yamir Moreno, and Silvio C. Ferreira.
    Efficient Gillespie algorithms for spreading phenomena in large and heterogeneous higher-order networks.
    Nature Communications (2026). DOI:10.1038/s41467-026-75402-0. https://www.nature.com/articles/s41467-026-75402-0
"""

print(banner)
