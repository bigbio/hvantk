# Installation

## Using Poetry (recommended)

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

## Using pip

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install -e .
```

## Prerequisites

- **Python** ≥ 3.10
- **Hail** ≥ 0.2.137
- **Java** 8 or 11 (required by Hail/Spark)
