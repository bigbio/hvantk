#!/bin/bash
#
# Setup fresh hvantk conda environment for the HGC scalability benchmark
#
# Usage:
#   bash setup_hvantk_env.sh [conda_env_name] [python_version]
#
# Default environment name: hvantk
# Default Python version: 3.10
#
# This script will:
# 1. Create a fresh conda environment (or recreate if exists)
# 2. Install all required dependencies
# 3. Install hvantk from the local repository
# 4. Verify everything works

set -euo pipefail

ENV_NAME="${1:-hvantk}"
PYTHON_VERSION="${2:-3.10}"

echo "========================================================================"
echo "Setting up fresh hvantk conda environment: $ENV_NAME"
echo "========================================================================"
echo ""

# Find hvantk repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HVANTK_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "hvantk repository: $HVANTK_ROOT"
echo "Python version: $PYTHON_VERSION"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found!"
    echo "Please install conda or make sure it's in your PATH"
    exit 1
fi

# Initialize conda for bash
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
else
    eval "$(conda shell.bash hook 2>/dev/null)" || true
fi

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "WARNING: Conda environment '$ENV_NAME' already exists!"
    echo ""
    read -p "Do you want to recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
        echo "✓ Existing environment removed"
        echo ""
    else
        echo "Using existing environment..."
        echo ""
    fi
fi

# Create the environment if it doesn't exist
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Step 1: Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y || {
        echo "ERROR: Failed to create conda environment"
        exit 1
    }
    echo "✓ Environment created"
    echo ""
else
    echo "Step 1: Using existing conda environment '$ENV_NAME'..."
    echo ""
fi

echo "Step 2: Activating environment..."
conda activate "$ENV_NAME"
echo "✓ Environment activated"
echo "  Python: $(which python)"
echo "  Python version: $(python --version)"
echo ""

echo "Step 3: Installing required packages..."
echo "  - hail (genomics analysis framework)"
echo "  - pandas (data manipulation)"
echo "  - numpy (numerical computing)"
echo "  - matplotlib (plotting)"
echo "  - seaborn (statistical visualization)"
echo ""

pip install --upgrade pip
pip install hail pandas numpy matplotlib seaborn || {
    echo "ERROR: Failed to install required packages"
    exit 1
}
echo "✓ Required packages installed"
echo ""

echo "Step 4: Installing hvantk from local repository..."
cd "$HVANTK_ROOT"

# Check if requirements.txt or pyproject.toml exists
if [ -f "requirements.txt" ]; then
    echo "  Found requirements.txt, installing dependencies..."
    pip install -r requirements.txt
fi

# Install hvantk in editable mode
pip install -e . || {
    echo "ERROR: Failed to install hvantk"
    exit 1
}
echo "✓ hvantk installed"
echo ""

echo "Step 5: Verifying installation..."
echo "  Checking hvantk imports..."
if python -c "import hvantk; print(f'  hvantk version: {hvantk.__version__ if hasattr(hvantk, \"__version__\") else \"unknown\"}')" 2>/dev/null; then
    echo "  ✓ hvantk module imports successfully"
else
    echo "  ✗ hvantk module failed to import"
    exit 1
fi

echo "  Checking hvantk.hgc functions..."
if python -c "from hvantk.hgc import combine_gvcfs, convert_vds_to_mt, compute_full_qc, convert_mt_to_multi_sample_vcf; print('  ✓ All HGC functions importable')" 2>/dev/null; then
    true
else
    echo "  ✗ hvantk.hgc functions failed to import"
    exit 1
fi

echo "  Checking hail..."
if python -c "import hail as hl; print(f'  Hail version: {hl.__version__}')" 2>/dev/null; then
    echo "  ✓ Hail imports successfully"
else
    echo "  ✗ Hail failed to import"
    exit 1
fi
echo ""

echo "========================================================================"
echo "✓ Setup complete!"
echo "========================================================================"
echo ""
echo "Fresh conda environment '$ENV_NAME' is ready for the HGC scalability benchmark."
echo ""
echo "Installed packages:"
conda list | grep -E "(hail|pandas|numpy|matplotlib|seaborn|hvantk)" || true
echo ""
echo "To run the benchmark:"
echo "  bash hgc_scalability_benchmark.sh --conda-env $ENV_NAME [other options]"
echo ""
echo "Or activate the environment and run:"
echo "  conda activate $ENV_NAME"
echo "  bash hgc_scalability_benchmark.sh [options]"
echo ""
echo "To deactivate when done:"
echo "  conda deactivate"
echo ""

