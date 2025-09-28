from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess
import sys
from pathlib import Path
import os

def run_fpm_command(cmd):
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        # Mostra stdout e stderr do comando que falhou
        raise RuntimeError(
            f"Command '{' '.join(cmd)}' failed with return code {e.returncode}\n"
            f"stdout:\n{e.stdout}\n"
            f"stderr:\n{e.stderr}"
        )

class InstallWithFPM(install):
    def run(self):
        # Check if the compiler is installed, via FC environment variable or default to gfortran
        compiler = os.environ.get("FC", "gfortran")
        fflags = os.environ.get("FFLAGS", "-O3 -ffree-line-length-512")

        try:
            subprocess.check_call([compiler, "--version"])
        except FileNotFoundError:
            sys.exit(f"❌ Fortran compiler '{compiler}' not found. Please install it or set the FC environment variable.")

        # Check if fpm is installed
        try:
            subprocess.check_call(["fpm", "--version"])
        except FileNotFoundError:
            sys.exit("❌ 'fpm' not found. Please install it before continuing: https://fpm.fortran-lang.org")

        # Compiles the Fortran executable
        print("🔧 Compiling executable with fpm...")

        # Directory where pip installs scripts/executables
        bin_dir = Path(sys.prefix) / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        # fpm install with check output and errors
        run_fpm_command(["fpm", "clean", "--all"])
        run_fpm_command([
            "fpm", "install",
            "--profile", "release",
            "--compiler", compiler,
            "--flag", fflags,
            "--prefix",
            str(Path(sys.prefix))
        ])

        # It will be installed
        print(f"✅ Executable installed at: {Path(sys.prefix) / 'bin' / 'hyperSIS_*'}")

        # Continues with the normal installation of the Python package
        super().run()

setup(
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=["fpm>=0.12.0"],
    cmdclass={
        "install": InstallWithFPM,
    },
)
