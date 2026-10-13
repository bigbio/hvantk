# GTEx fixture provenance and redistribution notice

`Liver.v11.eQTLs.signif_pairs.parquet` is a real 800-row slice of the GTEx v11
Liver cis-eQTL significant-pairs release (see `THIRD_PARTY_DATA.md` at the repo
root for the full source/licence entry).

GTEx's own data licence (https://gtexportal.org/home/license) requires anyone
redistributing GTEx data to carry a notice that the copy may not reflect the
most current/accurate release, and to acknowledge the source with a dated
statement. Any further redistribution of this fixture (or data derived from
it) should carry wording along these lines:

> This file is a small excerpt of the GTEx v11 Liver cis-eQTL release and may
> not reflect the most current or accurate data available from the GTEx
> Portal. Obtained from the GTEx Portal on 2026-10-13.

Checked 2026-10-13. GTEx's licence page is a JavaScript-rendered single-page
app that this pass could not fetch and re-quote byte-for-byte in this
environment; the requirement above is corroborated via independently-found
search-index snippets of the live page's acknowledgement-statement template
("The data used for the analyses described in this manuscript were obtained
from: … the GTEx Portal on MM/DD/YY and/or dbGaP accession number
phs000424.vN.pN on MM/DD/YYYY"), not a fully re-verified direct quote of every
clause. A maintainer with interactive browser access should confirm the exact
wording directly against gtexportal.org/home/license.
