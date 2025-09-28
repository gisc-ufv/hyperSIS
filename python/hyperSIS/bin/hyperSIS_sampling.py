import os
import sys
import subprocess
from pathlib import Path

def main_stub():
    bin_path = Path(__file__).parent / "hyperSIS_sampling"
    if not bin_path.exists():
        sys.exit("❌ Executable hyperSIS_sampling not found inside package!")
    os.execv(bin_path.as_posix(), [bin_path.as_posix()] + sys.argv[1:])
