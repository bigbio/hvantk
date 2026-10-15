# Reading the results

`pipeline_stats.json`:

| key | meaning | distrust the run when |
| --- | --- | --- |
| `n_query_samples`, `n_reference_samples` | sizes after loading | either is far from what you expected |
| `n_shared_variants` | variants present in both before filtering | below the minimum you accepted; the merge used a slice |
| `n_variants_after_filter`, `n_variants_after_ld_prune` | what the PCA saw | a few hundred or fewer at the continental level |
| `variance_explained_pc1`, `variance_explained_top5` | structure captured | the top components explain almost nothing: the shared variants carry no structure |
| `populations`, `n_populations` | the labels learned | a label with fewer than ten samples was present |
| `cv_accuracy` | cross-validated accuracy on the reference | 1.0 on a real panel (too few samples per label), or far below 0.9 at the continental level |
| `n_query_assigned`, `n_query_unassigned` | samples above and below `--min-prob` | most samples unassigned: wrong reference, admixture, or a build mismatch that survived the merge |

`predictions.tsv`: one row per sample, query first and the reference samples after
them (`_ancestry_source` says which); `predicted_ancestry` is the literal `unassigned`
for samples below the threshold, `ancestry_probability` the winning class's probability, `prob_<POP>` the full
vector, `_known_ancestry` for reference samples in the merged table, and the PCs. A
plot of PC1 against PC2 coloured by prediction, with the reference samples as the
background, is the one figure worth making; the HTML report contains it.

A run is reported with: the reference and its build, the shared-variant count, the
accuracy, the assigned fraction, and the threshold. A number without those is not
reproducible.
