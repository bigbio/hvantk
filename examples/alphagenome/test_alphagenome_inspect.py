#!/usr/bin/env python
"""Inspect the structure of a single AlphaGenome API response."""

import os
import sys

if not os.environ.get("ALPHAGENOME_API_KEY"):
    print("ERROR: ALPHAGENOME_API_KEY not set")
    sys.exit(1)

from alphagenome.data import genome as ag_genome
from alphagenome.models import dna_client as ag_client

api_key = os.environ["ALPHAGENOME_API_KEY"]
model = ag_client.create(api_key, timeout=120)

# Single variant: TP53 R175H
interval = ag_genome.Interval(chromosome="chr17", start=7150800, end=8199376)
variant = ag_genome.Variant(
    chromosome="chr17", position=7675088,
    reference_bases="C", alternate_bases="T",
)

print("Calling predict_variant for TP53 R175H...")
result = model.predict_variant(
    interval=interval,
    variant=variant,
    ontology_terms=["UBERON:0001157"],
    requested_outputs=[ag_client.OutputType.RNA_SEQ],
)

print(f"\nResult type: {type(result).__name__}")
print(f"Result dir: {[a for a in dir(result) if not a.startswith('_')]}")

if hasattr(result, "reference"):
    ref = result.reference
    print(f"\nreference type: {type(ref).__name__}")
    print(f"reference dir: {[a for a in dir(ref) if not a.startswith('_')]}")
    if hasattr(ref, "__dict__"):
        for k, v in vars(ref).items():
            if k.startswith("_"):
                continue
            print(f"  .{k}: type={type(v).__name__}", end="")
            if hasattr(v, "shape"):
                print(f", shape={v.shape}", end="")
            elif isinstance(v, (list, dict)):
                print(f", len={len(v)}", end="")
            elif isinstance(v, str):
                print(f", val={v[:80]!r}", end="")
            print()
            # If it's a dict or list, peek inside
            if isinstance(v, dict):
                for dk, dv in list(v.items())[:3]:
                    print(f"    [{dk!r}]: type={type(dv).__name__}", end="")
                    if hasattr(dv, "__dict__"):
                        print(f", attrs={[a for a in dir(dv) if not a.startswith('_')]}", end="")
                        for a in [a for a in dir(dv) if not a.startswith('_')][:5]:
                            av = getattr(dv, a)
                            if hasattr(av, "shape"):
                                print(f"\n      .{a}: shape={av.shape}, dtype={av.dtype}", end="")
                            elif callable(av):
                                pass
                            else:
                                print(f"\n      .{a}: type={type(av).__name__}", end="")
                    print()
            elif isinstance(v, list) and len(v) > 0:
                first = v[0]
                print(f"    [0]: type={type(first).__name__}", end="")
                if hasattr(first, "__dict__"):
                    print(f", attrs={[a for a in dir(first) if not a.startswith('_')]}", end="")
                print()

if hasattr(result, "alternate"):
    alt = result.alternate
    print(f"\nalternate type: {type(alt).__name__}")
    if hasattr(alt, "__dict__"):
        for k, v in vars(alt).items():
            if k.startswith("_"):
                continue
            print(f"  .{k}: type={type(v).__name__}", end="")
            if hasattr(v, "shape"):
                print(f", shape={v.shape}", end="")
            elif isinstance(v, (list, dict)):
                print(f", len={len(v)}", end="")
            print()

print("\nDone.")
