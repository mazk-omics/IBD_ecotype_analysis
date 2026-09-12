echo "===== executables ====="

find /home/mazekai/miniconda3/envs \
  -maxdepth 3 \
  \( -name python -o -name python3 -o -name R -o -name Rscript \) \
  -ls 2>/dev/null

echo
echo "===== omicverse ====="

/home/mazekai/miniconda3/bin/conda run -n omicverse python - <<'PY'
import sys
import importlib.metadata as md

print("Python:", sys.version)
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

echo
echo "===== R in environments ====="

for env in jupyterlab omicverse spatial; do
    echo "--- $env ---"
    if [ -x "/home/mazekai/miniconda3/envs/$env/bin/R" ]; then
        "/home/mazekai/miniconda3/envs/$env/bin/R" --version | head -1
    else
        echo "R: not installed"
    fi
done