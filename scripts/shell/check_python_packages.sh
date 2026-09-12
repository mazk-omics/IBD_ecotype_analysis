for env in jupyterlab omicverse spatial; do
    echo
    echo "========================================"
    echo "ENV: $env"
    echo "========================================"

    PY="/home/mazekai/miniconda3/envs/$env/bin/python"

    "$PY" - <<'PY'
import sys
import importlib.metadata as md

print("Python:", sys.version.split()[0])
print("Executable:", sys.executable)

packages = [
    "numpy",
    "pandas",
    "scipy",
    "pyarrow",
    "anndata",
    "scanpy",
    "h5py",
    "matplotlib",
]

for p in packages:
    try:
        print(f"{p:12s} {md.version(p)}")
    except Exception:
        print(f"{p:12s} NOT INSTALLED")
PY

done