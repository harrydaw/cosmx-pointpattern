"""Execute a notebook in-place with the repo venv kernel (bypasses the stale
global kernelspec). Usage: python docs/_run_nb.py notebooks/13_diagnostics.ipynb
"""
import sys, os
from pathlib import Path
import nbformat
from nbclient import NotebookClient
from jupyter_client.manager import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

nb_path = Path(sys.argv[1]).resolve()
nb = nbformat.read(str(nb_path), as_version=4)

# KernelManager whose spec points at THIS interpreter
km = KernelManager()
km.kernel_name = "python3"
spec = km.kernel_spec  # may raise if missing; override argv regardless
spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]

client = NotebookClient(
    nb,
    timeout=1800,
    kernel_manager=km,
    resources={"metadata": {"path": str(nb_path.parent)}},
)
client.execute()
nbformat.write(nb, str(nb_path))
print("EXECUTED", nb_path.name)
