#!/usr/bin/env python3
"""Cut a tiny real-data fixture from an AlphaGenome tidy-scores shard, for the
hvantk `alphagenome:predictions` plugin test suite.

This script records exactly how `clinvar-subset.parquet` was produced; it is not
a general-purpose tool. AlphaGenome's licence requires the sibling NOTICE.md to
travel with the output. Re-running it on the same source shard reproduces the
subset byte for byte.

Input
-----
The source shard: AlphaGenome scores for 411 ClinVar canonical-splice variants
from an AlphaGenome SDK `score_variant` -> `variant_scorers.tidy_scores()` run on
2026-06-24 (RECOMMENDED scorers, 1 Mb interval, human, API backend), written
with pandas/pyarrow. 123,492,790 bytes; 16,642,316 rows.
sha256 (verified 2026-10-06):
261cf87273f1b36d20db640305f414b2e9bdc55824cddff9e4cc4bcc4c6745a1

Variants
--------
The three ClinVar variants of the first cut of this fixture, pinned here instead
of re-drawn so the fixture keeps its keys:

- chr3:39408741:T>C  -- has SPLICE_JUNCTIONS rows
- chr6:112216367:C>A -- negative raw_score on a gene-level scorer
- chrX:153694448:T>G -- on chrX

Rows
----
For every (variant, output_type, variant_scorer) group -- one group per summary
field the builder emits -- the subset keeps:

1. the row with the largest |quantile_score|, ties broken by the builder's own
   rule (track_name, then gene_id, then quantile_score, then raw_score,
   ascending, missing last), so it is the row the builder's `top_*` fields
   report from the full shard;
2. the row with the largest |raw_score|, so `max_abs_raw` computed from the
   fixture equals `max_abs_raw` computed from the full shard;
3. one random row (SEED).

Per variant it also keeps every heart SPLICE_SITE_USAGE row: gtex_tissue
containing "Heart", or ontology_curie in HEART_CURIES. The curies also tag an
ENCODE track, so the two sets differ (SKILL.md s 4). With both kept, a heart
build of the fixture (`ontology_curies` = HEART_CURIES) equals the heart build
of the full shard, which the round-trip test pins.

Every SPLICE_JUNCTIONS row in the shard carries junction_Start/junction_End, so
the SPLICE_JUNCTIONS groups supply rows with the junction columns set; this is
asserted below.

Row order
---------
The subset is written in reverse source order. In source order the builder's
pick came first in every tie at the top |quantile_score| that the subset holds,
so a builder that ignored the tie-break and kept the first row it read still
matched the snapshot. Reversed, it does not: the script asserts that no tied
group starts with the builder's pick.

Seed
----
SEED = 20261006, fixed. Used only to draw the random row of step 3.

Run
---
Needs numpy, pandas and pyarrow:

    python3 make_fixture.py <source-shard.parquet> [--out <subset.parquet>]

Output
------
Writes the subset (by default clinvar-subset.parquet next to this script) and
prints its row counts per variant and per scorer, the tied groups checked, its
row count and file size, and the sha256 of both the source shard and the subset.
"""

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SEED = 20261006
OUT = Path(__file__).resolve().parent / "clinvar-subset.parquet"

VARIANTS = ("chr3:39408741:T>C", "chr6:112216367:C>A", "chrX:153694448:T>G")
HEART_CURIES = ("UBERON:0006566", "UBERON:0006631")
GROUP = ["variant_id", "output_type", "variant_scorer"]
# The builder's top-row order after |quantile_score|; see the module docstring.
TIE_BREAK = ["track_name", "gene_id", "quantile_score", "raw_score"]


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ties_led_by_the_pick(df):
    """Count the tied groups, and list those whose first row (in file order) at the
    top |quantile_score| is also the builder's pick."""
    d = df.assign(pos=np.arange(len(df)), abs_q=df["quantile_score"].abs())
    top = d[d["abs_q"] == d.groupby(GROUP)["abs_q"].transform("max")]
    tied = top[top.groupby(GROUP)["pos"].transform("size") > 1]
    picks = tied.sort_values(TIE_BREAK, na_position="last", kind="mergesort")
    picks = picks.groupby(GROUP).head(1)
    firsts = tied.sort_values("pos").groupby(GROUP).head(1)
    led = picks.merge(firsts[GROUP + ["pos"]], on=GROUP + ["pos"])
    return tied.groupby(GROUP).ngroups, led[GROUP].values.tolist()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("shard", help="source tidy_scores() parquet (read-only)")
    parser.add_argument("--out", type=Path, default=OUT, help="subset to write")
    args = parser.parse_args()

    rng = np.random.default_rng(SEED)

    print(f"source shard sha256: {sha256_of(args.shard)}")
    table = pq.read_table(args.shard)
    print(f"source rows: {table.num_rows:,}")

    cols = table.select(
        GROUP + TIE_BREAK + ["gtex_tissue", "ontology_curie"]
    ).to_pandas()
    cols["row_idx"] = np.arange(table.num_rows)
    sub = cols[cols["variant_id"].isin(VARIANTS)]
    absent = sorted(set(VARIANTS) - set(sub["variant_id"]))
    assert not absent, f"variants not in the shard: {absent}"

    keep: set[int] = set()

    # 1. The builder's top row per group.
    ranked = sub.assign(neg_abs_q=-sub["quantile_score"].abs()).sort_values(
        GROUP + ["neg_abs_q"] + TIE_BREAK, na_position="last", kind="mergesort"
    )
    keep.update(ranked.groupby(GROUP, sort=True).head(1)["row_idx"].tolist())

    # 2. The largest |raw_score| per group.
    abs_raw = sub["raw_score"].abs()
    top_raw = abs_raw.groupby([sub[c] for c in GROUP], sort=True).idxmax()
    keep.update(sub.loc[top_raw.values, "row_idx"].tolist())

    # 3. One random row per group, groups visited in sorted order.
    for _, group in sub.groupby(GROUP, sort=True):
        keep.add(int(group["row_idx"].iloc[rng.integers(len(group))]))

    # Every heart SPLICE_SITE_USAGE row, under either heart definition.
    ssu = sub[sub["output_type"] == "SPLICE_SITE_USAGE"]
    by_tissue = ssu["gtex_tissue"].fillna("").str.contains("Heart", case=False)
    heart = by_tissue | ssu["ontology_curie"].isin(HEART_CURIES)
    keep.update(ssu.loc[heart, "row_idx"].tolist())
    print(f"heart SPLICE_SITE_USAGE rows kept: {int(heart.sum())}")

    # Reverse source order, so that no tie is won by the first row read.
    order = sorted(keep, reverse=True)
    subset = table.take(pa.array(order, type=pa.int64()))
    print(f"\nsubset rows: {subset.num_rows}")

    source_schema = table.schema.remove_metadata()
    subset_schema = subset.schema.remove_metadata()
    assert subset_schema.equals(source_schema), (
        "subset schema diverged from source schema"
    )
    print(
        "schema check: subset.schema.remove_metadata() == source.schema.remove_metadata() -> OK"
    )

    df = subset.to_pandas()
    sj = df[df["output_type"] == "SPLICE_JUNCTIONS"]
    assert sj["junction_Start"].notna().any() and sj["junction_End"].notna().any(), (
        "no SPLICE_JUNCTIONS row with junction columns set"
    )
    print(
        f"SPLICE_JUNCTIONS rows: {len(sj)}, with junction_Start/End set: "
        f"{int((sj['junction_Start'].notna() & sj['junction_End'].notna()).sum())}"
    )
    n_tied, led = ties_led_by_the_pick(df)
    assert not led, f"tied groups that start with the builder's pick: {led}"
    print(f"tied groups: {n_tied}, none starts with the builder's pick")
    print("\nrows per variant:")
    print(df.groupby("variant_id").size().to_string())
    print("\ndistinct (output_type, variant_scorer) groups per variant:")
    print(df.groupby("variant_id")["variant_scorer"].nunique().to_string())
    print("\nrows per (output_type, variant_scorer):")
    print(df.groupby(["output_type", "variant_scorer"]).size().to_string())

    pq.write_table(subset, args.out)
    print(f"\nwrote {args.out}: {args.out.stat().st_size:,} bytes")
    print(f"subset sha256: {sha256_of(args.out)}")


if __name__ == "__main__":
    main()
