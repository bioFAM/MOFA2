# Enrichment analysis on factor weights

Gene set enrichment on MOFA weights asks: for a given factor, are the
high-weighted features concentrated in any known gene set? It converts a ranked
feature list into a biological statement, and it is usually the step that turns
"Factor 3 exists" into "Factor 3 is an interferon response".

## Feature sets

`run_enrichment()` requires a **binary matrix**, gene sets in rows, features in
columns, 1 where the feature belongs to the set. This is stricter than most GSEA
tooling — a list of character vectors will be rejected.

The GSEA tutorial uses the pre-processed sets shipped in `MOFAdata`:

```r
library(MOFAdata)
data("reactomeGS")             # Reactome (human)
data("MSigDB_v6.0_C2_human")   # curated pathway sets (C2)
data("MSigDB_v6.0_C5_human")   # Gene Ontology (C5); _mouse variants also exist
feature.sets <- MSigDB_v6.0_C5_human
```

The 10x Multiome tutorial instead queries MSigDB live with `msigdbr` and reshapes
it to the same binary matrix (one row per set, one column per gene) — useful when
you want a category or species the frozen `MOFAdata` sets don't cover.

Converting your own sets:

```r
genes <- unique(unlist(gene_set_list))
feature.sets <- t(sapply(gene_set_list, function(g) as.integer(genes %in% g)))
colnames(feature.sets) <- genes
```

**Feature name matching is the usual failure point.** Column names of
`feature.sets` must match the feature names in the view exactly — same
identifier type (symbol vs. Ensembl), same case. MOFA also appends a view suffix
to features present in several views, which must be stripped first:

```r
colnames(feature.sets) <- toupper(colnames(feature.sets))
features_names(model)[["mRNA"]] <- toupper(features_names(model)[["mRNA"]])
```

Check the overlap before running — a near-empty intersection produces results
that look valid and mean nothing:

```r
length(intersect(colnames(feature.sets), features_names(model)[["mRNA"]]))
```

## Running it

```r
res.positive <- run_enrichment(model, view = "mRNA", factors = "all",
                               feature.sets    = feature.sets,
                               sign            = "positive",
                               statistical.test = "parametric",
                               set.statistic   = "mean.diff",
                               min.size        = 10,
                               p.adj.method    = "BH", alpha = 0.1)

res.negative <- run_enrichment(model, view = "mRNA", factors = "all",
                               feature.sets = feature.sets, sign = "negative")
```

**Run the two signs separately.** A factor's positive and negative weights
describe opposite biology, and merging them (taking absolute values) can dilute
the signal when the two directions belong to different pathways. Both the GSEA and
10x tutorials recommend running `sign = "positive"` and `sign = "negative"`
separately — optionally also jointly (`sign = "all"`) — and reading each against
the direction of the factor.

### Choosing the test

`statistical.test`:

- `"parametric"` — a t-test on the weights; effectively instant. What both the
  GSEA and 10x tutorials use by default.
- `"cor.adj.parametric"` — adjusts for correlation between features. The GSEA
  tutorial shows it gives **more conservative** p-values than the plain parametric
  test, though the two are highly correlated; it costs ~10 min where parametric is
  near-instant.
- `"permutation"` — non-parametric and the most conservative, and much the
  slowest (`nperm = 1000` default). Reasonable to confirm a headline result.

The tutorial's own comparison is parametric vs. cor.adj.parametric, and its
takeaway is that the adjustment shifts p-values conservatively without
reordering them much — so parametric is a fine default for exploration, with the
adjusted test as a check.

`set.statistic` is `"mean.diff"` (default) or `"rank.sum"`. `"mean.diff"` uses
weight magnitudes; `"rank.sum"` only their ordering, so it is more robust to a few
extreme weights.

## Visualising

```r
plot_enrichment(res.positive, factor = 1, max.pathways = 15)

plot_enrichment_heatmap(res.positive, alpha = 0.1)   # all factors at once

plot_enrichment_detailed(res.positive, factor = 1,
                         max.genes = 8, max.pathways = 5)
```

`plot_enrichment_heatmap()` is the best overview — it shows which factors have
coherent biology at all. `plot_enrichment_detailed()` shows which individual
genes drive an enriched set, which is the check on whether a hit rests on one or
two genes.

## Motif / TF enrichment for ATAC views

For chromatin views the question becomes: which transcription-factor motifs are
enriched among the peaks driving a factor? The 10x Multiome tutorial does this
**with the same `run_enrichment()` function** — the only change is swapping the
gene-set matrix for a **motif × peak** binary matrix (1 where a peak contains a
motif). There is no separate motif-enrichment step; motif enrichment is GSEA with
a different annotation matrix.

The motif–peak matrix comes from the Signac motif object attached to the ATAC
assay (built from JASPAR2020 PWMs via `CreateMotifMatrix` upstream). Signac stores
it as peaks × motifs, so transpose it into the sets × features orientation
`run_enrichment()` expects:

```r
# motifs x peaks, i.e. feature-sets x features
motif.matrix <- t(as.matrix(seurat[["ATAC_distal"]]@motifs@data))

motif.enrichment.positive <- run_enrichment(model,
  view = "ATAC_distal", factors = 1:2,
  feature.sets = motif.matrix, sign = "positive")

motif.enrichment.negative <- run_enrichment(model,
  view = "ATAC_distal", factors = 1:2,
  feature.sets = motif.matrix, sign = "negative")
```

Three practical notes from the tutorial: run it on the ATAC view carrying the most
factor variance (there, `ATAC_distal` — distal peaks explained more than
promoters); keep the +/- sign separation, since peaks opening and closing along a
factor implicate different TFs; and peak-name matching bites here as it does for
genes — the tutorial rewrites `:` to `-` in the feature names so they line up with
the motif-matrix columns.

## Interpreting results honestly

Enrichment is a hypothesis generator. Two caveats worth passing on: gene sets
overlap heavily, so several "independent" hits often reflect one biological
program; and a factor with no enrichment is not meaningless — it may be driven by
biology with no annotated set, which is often the more interesting case.