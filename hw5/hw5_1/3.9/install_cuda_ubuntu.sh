#!/usr/bin/env bash
set -euo pipefail

echo "[info] Check OS"
source /etc/os-release
echo "Ubuntu version: ${VERSION_ID:-unknown}"

echo "[info] Check NVIDIA GPU"
if ! lspci | grep -i nvidia >/dev/null 2>&1; then
  echo "[warn] No NVIDIA GPU detected by lspci."
  echo "[warn] CUDA requires a CUDA-capable NVIDIA GPU."
  exit 1
fi

echo "[info] Check current driver"
nvidia-smi || true

echo "[info] Install basic build tools"
sudo apt-get update
sudo apt-get install -y build-essential gcc g++ make dkms linux-headers-$(uname -r)

echo
echo "[next step]"
echo "Open the official NVIDIA CUDA Linux installation guide:"
echo "  https://docs.nvidia.com/cuda/cuda-installation-guide-linux/"
echo
echo "Choose your exact Ubuntu version and preferred method:"
echo "  - deb(network)"
echo "  - deb(local)"
echo "  - runfile"
echo
echo "[common verification after install]"
echo "  nvcc --version"
echo "  nvidia-smi"
echo
echo "[build GPU version after CUDA install]"
echo "  make fractal_gpu"