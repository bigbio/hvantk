# NOTICE: dbNSFP data

**`dbNSFP4_v49a_example_variants.bgz` is an excerpt of dbNSFP v4.9a (academic branch),
licensed under CC BY-NC-ND 4.0 (<https://creativecommons.org/licenses/by-nc-nd/4.0/>;
terms: <https://www.dbnsfp.org/license/>). Non-commercial use only, with attribution;
do not share modified versions. These data are not covered by this repository's MIT
licence.**

Commercial use requires dbNSFP's paid commercial licence. The CADD, VEST, M-CAP,
MutScore, PolyPhen-2, PrimateAI and RGC Million Exome scores in the academic branch
also need commercial licences from their authors.

## What this is

The header (all 458 columns) and the first 4,999 variants on chromosome 10 of the
dbNSFP v4.9a variant table. Once decompressed it is byte-identical to the release;
nothing was edited, trimmed or reordered. It is stored BGZF-compressed.

## Attribution

dbNSFP, © 2024–2026 Genos Bioinformatics LLC, <https://www.dbnsfp.org>. Cite: Liu X, Li C,
Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific
functional predictions and annotations for human nonsynonymous and splice-site SNVs.
*Genome Medicine* 12, 103 (2020). <https://doi.org/10.1186/s13073-020-00803-9>

## Regenerating

Extract whole lines unedited (`hvantk/skills/dbnsfp/SKILL.md` § 8). Trimming columns
or editing values would make a modified version, which the licence does not allow
sharing. Record any change in the dbNSFP entry of `THIRD_PARTY_DATA.md`.
