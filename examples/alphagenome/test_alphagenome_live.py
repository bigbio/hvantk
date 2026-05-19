#!/usr/bin/env python
"""Live smoke test for AlphaGenome streamer with real API calls."""

import json
import logging
import os
import sys
import tempfile

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VARIANTS_TSV = os.path.join(SCRIPT_DIR, "test_alphagenome_variants.tsv")
CONFIG_YAML = os.path.join(SCRIPT_DIR, "alphagenome_test_config.yaml")

# Verify API key is set
if not os.environ.get("ALPHAGENOME_API_KEY"):
    print("ERROR: ALPHAGENOME_API_KEY environment variable not set")
    sys.exit(1)

from hvantk.core.streamers.alphagenome import (
    AlphaGenomeStreamer,
    load_config,
    _load_variants_from_tsv,
    compute_intervals,
)

# --- Step 1: Validate config loading (picks up env var) ---
print("\n=== Step 1: Config loading ===")
config = load_config(CONFIG_YAML)
print(f"  API key resolved: {'*' * 8}{config['api']['key'][-4:]}")
print(f"  Ontology terms: {config['ontology']['terms']}")
print(f"  Output types: {config['ontology']['output_types']}")
print(f"  Adaptive intervals: {config['intervals']['adaptive']}")

# --- Step 2: Validate variant loading ---
print("\n=== Step 2: Variant loading ===")
variants = _load_variants_from_tsv(VARIANTS_TSV)
for v in variants:
    print(f"  {v.chrom}:{v.pos} {v.ref}>{v.alt}")

# --- Step 3: Validate interval computation ---
print("\n=== Step 3: Interval computation ===")
intervals = compute_intervals(variants, config)
print(f"  {len(variants)} variants grouped into {len(intervals)} intervals:")
for iv, vs in intervals:
    print(f"    {iv.chrom}:{iv.start}-{iv.end} ({len(vs)} variant(s))")

# --- Step 4: Run full streamer with real API ---
print("\n=== Step 4: Running AlphaGenome streamer (live API) ===")
with tempfile.TemporaryDirectory(prefix="alphagenome_test_") as output_dir:
    streamer = AlphaGenomeStreamer(
        input_path=VARIANTS_TSV,
        output_dir=output_dir,
        config_path=CONFIG_YAML,
        no_resume=True,
        chunk_size=10,
    )

    streamer.setup()

    all_results = {}
    for batch in streamer.stream():
        for key, result in batch.items():
            all_results[key] = result
            print(f"  Prediction received: {key}")

    streamer.teardown()

    # --- Step 5: Check outputs ---
    print(f"\n=== Step 5: Results summary ===")
    print(f"  Total predictions: {len(all_results)}")

    predictions_path = os.path.join(output_dir, "predictions.json")
    if os.path.exists(predictions_path):
        with open(predictions_path) as f:
            saved = json.load(f)
        print(f"  Saved to predictions.json: {len(saved)} entries")

        # Show structure of first prediction
        first_key = next(iter(saved))
        first_val = saved[first_key]
        print(f"\n  Sample prediction ({first_key}):")
        print(f"    Keys: {list(first_val.keys())}")
        for k, v in first_val.items():
            if isinstance(v, dict):
                print(f"    {k} sub-keys: {list(v.keys())}")
            elif isinstance(v, str):
                print(f"    {k}: {v[:200]}...")
    else:
        print("  WARNING: predictions.json not found!")

    checkpoint_dir = os.path.join(output_dir, "_checkpoints")
    if os.path.isdir(checkpoint_dir):
        state_path = os.path.join(checkpoint_dir, "state.json")
        if os.path.exists(state_path):
            with open(state_path) as f:
                state = json.load(f)
            print(f"\n  Checkpoint state:")
            print(f"    Completed intervals: {len(state['completed_intervals'])}")
            print(f"    Failed variants: {len(state['failed_variants'])}")
            if state['failed_variants']:
                for fv in state['failed_variants']:
                    print(f"      FAILED: {fv['chrom']}:{fv['pos']} {fv['ref']}>{fv['alt']} - {fv['reason']}")

print("\n=== Done ===")
