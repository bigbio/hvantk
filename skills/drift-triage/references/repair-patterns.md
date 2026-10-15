# Repair patterns for a schema-signal drift

One pattern per shape of change, in the order a triage usually meets them. Each
plugin's rename map lives at `hvantk/skills/<provider>/shared/constants.py`
(for example `CLINGEN_GENE_DISEASE_FIELDS`, `HGNC_GENE_FIELDS`); a builder imports
EVERY column in the raw file and renames only the ones the map declares --
`hvantk/skills/_conventions/SKILL.md` section 5 has the builder contract, read it
there, it is not restated here.

## Renamed column

The new file has a column the rename map does not know, and is missing the old name.
Add the new source name to the map alongside the old one -- keep both rather than
replace, so a build from an older archived raw file still resolves:

```python
CLINGEN_GENE_DISEASE_FIELDS = {
    "GENE ID (HGNC)": "hgnc_id",   # kept for older raw files
    "HGNC ID": "hgnc_id",          # the new name
    ...
}
```

If the renamed column is the dataset's key (clingen's `hgnc_id`, keyed by
`(hgnc_id, mondo_id)`), an unrepaired build fails loudly at `.key_by(...)` with the field
not found -- that failure is itself the signal that this column is the key, not a
cosmetic one.

## Removed column

- **Feeds no output field** (nothing in the artifact's schema is built from it):
  nothing to repair. A builder's existing `{k: v for k, v in FIELDS.items() if k in
  get_row_fields(ht)}` filtering already drops a rename entry whose source column is
  gone, and the column disappears from the output along with it.
- **Feeds an output field**: stop and decide with the user before touching code. The
  output contract is changing -- the field either disappears from every downstream join
  or gets a documented fallback -- and that is a product decision, not a mechanical fix
  this skill makes on its own.

## Added column

A column the rename map does not declare is NOT dropped: the builder's `.rename(...)`
only renames the keys it is given and leaves every other column exactly as it is, so an
added upstream column survives in the output under its own, un-renamed, upstream header
text. This is why a committed schema snapshot -- built once from a small fixture -- is a
FLOOR on a full build's schema, not an exact description of it: upstream can and does
add columns the fixture never saw, and they show up anyway. Nothing to repair unless the
new column should be mapped to a proper field name, which is a product decision like
the removed-feeds-an-output-field case above, just in the opposite direction.

## Type or format change in a column

The rename map only renames; it does not coerce. Most builders import with
`impute=False`, so every field stays a string regardless of what it looks like upstream
-- a format change (a date format, a delimiter inside a pipe-separated field) usually
shows up as a changed VALUE when comparing a fresh build against your own previous one,
not as a schema diff the drift probe's header-only fingerprint can see. Diff the new
file against the fixture by eye for this one.

## Header marker or metadata lines moved

ClinGen's builder finds the header by scanning for the first line containing
`GENE SYMBOL` (`hvantk/skills/clingen/builder.py`); if upstream moves the metadata block
or renames that marker, the builder's own `RuntimeError: ClinGen header "GENE SYMBOL"
not found` is the signal. Fix the marker string the builder scans for, not the rename
map -- this is a preprocessing change, not a column rename.

## Fixture refresh

Only when the committed fixture is no longer representative (a tested value was
withdrawn, a new column needs exercising). Sub-sample a representative set of rows from
the live file and regenerate both the fixture and its schema/row snapshots together --
do NOT assume the keys the existing `SAMPLE_KEYS` names are still a reasonable target:
those are specific curated entries, and a real build is not guaranteed to still carry
them (a real ClinGen build has been seen to raise `KeyError` for a key the committed
fixture's own `SAMPLE_KEYS` lists, because that exact gene-disease pairing was not in
the live data at the time). Pick new keys from whatever the live file actually has. State the refresh explicitly in the commit message -- a silent fixture refresh and
a silent column drop both leave only the JSON diff for a reviewer to notice, and look
the same from there.

## When snapshots may be regenerated

`pytest hvantk/skills/<provider>/tests -m hail --regenerate-snapshots`, or
`templates/run-plugin-tests.py --repo "$HVANTK_REPO" --dataset <provider:dataset> -- --regenerate-snapshots`,
only after the repair, and only to reflect an INTENTIONAL change (a new optional column
entering the schema, a widened array, a refreshed date value) -- never to make a build
match by editing the JSON by hand, which hides whatever the contract actually changed
to. State what changed and why in the commit; a reviewer reads the schema and row diff
against that statement, not the other way around.
