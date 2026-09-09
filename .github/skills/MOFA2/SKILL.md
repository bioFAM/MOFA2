---
name: mofa2-analysis
description: >-
  Guides multi-omics factor analysis with the MOFA2 R package (MOFA and MEFISTO):
  judging whether MOFA fits the data at all, profiling the data to propose sensible
  model configurations, checking that each omics layer is ready for the model,
  training it, and interpreting factors, weights and variance decomposition. Use this
  whenever the user mentions MOFA, MOFA2, MEFISTO, mofapy2, multi-omics integration,
  or latent factor analysis across several assays — and also when they describe
  having several molecular data types measured on the same samples (RNA + methylation
  + drug response, scRNA + scATAC, CITE-seq, multi-modal microbiome, spatial or
  time-course omics) and want to find shared axes of variation, even if they never
  say the word "MOFA".
---

# MOFA2 analysis

MOFA is an unsupervised factor analysis model for data where several molecular
assays were measured on the same samples. It decomposes that data into a small
number of latent factors, each carrying a weight for every feature in every
assay. The output is an ordination — like PCA, but sharing one latent space
across modalities, and reporting how much variance each factor explains in each
assay.

Your job is to help someone get a *trustworthy* MOFA model and then read it
correctly.

## Where your remit ends

You advise on **MOFA**: whether it fits, and how to configure, build, train and
read the model. How the data is normalised, transformed, batch-corrected or
filtered on the way in is the analyst's call — it depends on the assay, the
platform and upstream processing, which they know and a matrix doesn't reveal.
Hold that line even when a transform looks obvious; your suggestions get adopted,
so a wrong one is expensive.

Measure instead. Run the diagnostics, report what they found, name the
consequence for the model ("a leading axis tracking per-sample totals will come
back as Factor 1"), ask what state each view is in — then let the user decide. If
they ask you to implement a step they have chosen, do it without reopening the
choice.

## Two postures, and why they differ

**Before training: be cautious and consultative.** Almost every bad MOFA result
traces back to a decision made before `run_mofa` was ever called — data that
wasn't normalised the way this analysis needed, a batch effect left in, one view
swamping the others, groups defined for the wrong reason. These mistakes don't
throw errors. They produce a model that trains cleanly, plots beautifully, and
answers a different question than the user asked. By the time that is visible,
hours of compute and days of interpretation have gone into it.

So don't write training code until you know what you're working with. Profile the
data (below), establish the facts in the gate, ask about what you genuinely can't
determine, and say what you're assuming when you proceed anyway.

**After training: be proactive.** Once the model exists it is fixed, and
inspecting it is cheap and reversible. Here, offering the next diagnostic
unprompted is a service rather than a risk. Follow the characterisation loop,
suggest what to look at next, and interpret what comes back.

Throughout, **generate code for the user to run rather than running it
yourself**, unless they ask otherwise. Their data is usually large, often on a
remote machine, and what they need is a script they can rerun, edit and keep.

## Stage 0 — Is MOFA right, and which variant?

Check this first. Steering someone away from MOFA when it doesn't fit is more
valuable than helping them run it well.

MOFA needs **~15+ samples minimum** — it is a variance-decomposition method, and
below that there is no variance structure to decompose — and **substantially
matched samples across assays**. Missing values are fine: MOFA drops them from
the likelihood with no hidden imputation step, and matrix factorisation is
robust to a lot of them. Entirely disjoint sample sets are not fine.

Then pick the variant:

| Signal in the user's data or description | Variant |
|---|---|
| Several assays, no ordering among samples | **plain MOFA** |
| Samples carry a *continuous* covariate: time, developmental stage, spatial coordinates | **MEFISTO** → `references/mefisto.md` |
| Samples fall into predefined categories, and the question is *which sources of variation are shared vs. group-specific* | **multi-group MOFA** |
| Samples fall into categories and the question is *what separates the groups* | **NOT multi-group** — see below |

The multi-group distinction is the most misunderstood part of MOFA, so state it
plainly whenever groups come up. The multi-group framework centers features
*within* each group before fitting. It therefore cannot, by construction, find a
factor that separates the groups — that information has been removed. Its
purpose is to ask whether the same axes of variation operate in each group. A
user who wants group separation should run ordinary MOFA on all samples and test
factors against the group label afterwards; if a factor separates the groups,
MOFA will find it unsupervised.

The same logic covers known covariates generally. Don't put age, sex or treatment
into the model as covariates. Learn factors unsupervised, then relate them to
covariates with `correlate_factors_with_covariates`. Covariates supplied to the
model tend to be dropped in favour of a latent version of themselves anyway,
because the measured label rarely matches the molecular reality.

## Stage 1 — Profile the data and propose configurations

This is what makes the skill more useful than the documentation: rather than
asking the user to specify everything, look at what they have and *propose*.
Data shape determines several settings almost deterministically, so derive them
instead of interrogating.

Run the bundled profiler, which reports dimensions, missingness, view balance,
candidate MEFISTO covariates, likelihood hints and GPU availability, and prints
a suggested configuration with the reasoning attached:

```bash
Rscript scripts/profile_mofa_data.R --help
```

It accepts a list of matrices saved as `.rds`, a trained/untrained MOFA object,
or an existing `.hdf5` model, plus an optional sample-metadata table. If running
it isn't practical (data on a cluster, sensitive data), apply the same logic by
hand from the table in `references/options.md` — the profiler encodes exactly
those rules.

The signals it acts on, and what each implies:

| Observation | Proposal | Why |
|---|---|---|
| N < 15 | Warn: MOFA likely inappropriate | Not enough samples to decompose variance |
| N ≤ 25 / ≤1e3 / ≤1e4 / >1e4 | `num_factors` 5 / 15 / 20 / 25 | The package's own sample-size heuristic |
| N > ~1e4 | Consider `stochastic = TRUE` | Full-batch updates get expensive; SVI subsamples |
| N > 1e5 | `use_float32 = TRUE` | Halves memory; MOFA2 does this automatically |
| Widest view ≫ narrowest (>10×) | Report the ratio; `weight_views = TRUE` is the in-model lever | Wide views dominate the factors and crowd out small ones; reducing the wide view instead is the analyst's call |
| Numeric metadata column with many distinct values, especially named time/age/stage/day/x/y | Candidate MEFISTO covariate — **offer, don't assume** | Continuous structure is what MEFISTO exploits |
| MEFISTO with many distinct covariate values (>~1000) | `sparseGP = TRUE`, tune `frac_inducing` | Exact GP inference scales cubically in covariate points |
| `cupy` importable and an NVIDIA GPU present | Offer `gpu_mode = TRUE` | Large models train far faster; only worth mentioning if actually available |
| A view is non-negative integers (counts) | Report it; ask whether counts are the intended model input | Gaussian on raw counts is a poor fit, and poisson is rarely used on real data — but which input to pass is the user's decision |
| A view is narrow and binary (e.g. mutations) | `bernoulli` | Directly appropriate — this is what the CLL tutorial does; sparsity is fine |
| A view is wide and binary/near-binary (e.g. scATAC peaks) | `bernoulli` works but is slow at this width — say so and let the user weigh it | The cost is a property of the model; whether to supply a continuous version of the view instead is theirs |
| Features with no observation in any group | Recommend dropping | They contribute nothing and trigger warnings |

The three likelihood rows in that table turn on measured width and sparsity rather
than on the modality's name — two views that both "look binary" get opposite
answers depending on how wide they are, so run `recommend_likelihoods.R` instead
of eyeballing it. Likelihood is one of the few pre-training choices genuinely
yours to advise on, because it is a model option rather than a data
transformation. Note that automatic likelihood inference is currently disabled in
the package, so every non-gaussian view must be set explicitly.

Present proposals as proposals. Say what you observed, what you'd set, and why —
then let the user override. A MEFISTO suggestion in particular changes the whole
analysis, so never silently switch them onto it because a column looked
continuous.

## Stage 2 — The pre-training gate

The profiler settles data *shape*. It cannot settle data *provenance*, and that
is where the expensive mistakes live. Establish these, inferring what you can
from what the user has said, from their code, or from the data itself — then ask
in **one batched round**, not a drip of questions, about whatever is still
undetermined *and* would change your recommendation.

Items 1 and 2 are the ones most worth asking about: they cause the most damage
and users volunteer them least.

1. **What each view is, and its current normalisation state.** "Counts",
   "log-CPM", "VST", "beta values", "CLR-transformed" behave completely
   differently in the model, and "RNA-seq data" alone does not distinguish them.
   Record it and set the likelihood accordingly; it is context you need, not an
   opening to propose a different transform. → `references/data-prep.md`
2. **Known batch or technical covariates, and whether they've been handled.**
   MOFA has no batch term, so anything left in the data can become a factor —
   say so plainly, since users often assume the model accounts for it. If they do
   intend to remove technical structure, flag the ordering: it has to happen
   *before* feature selection, or the batch effect distorts which features look
   variable and gets baked into the feature set.
3. **Whether groups are intended, and for what question** → apply the Stage 0
   test.
4. **What the user actually wants from the model** — major axes of biological
   variation, imputation, features for a downstream classifier? This changes
   `num_factors` more than any data property does.

When you must proceed without an answer, take the conventional default, state it
in one line, and flag which downstream conclusions depend on it. Don't stall a
whole analysis on one unknown.

## Stage 3 — Check readiness and report

The highest-leverage stage, and the one where your role is narrowest: you
measure, the analyst decides. `scripts/check_preprocessing.R` flags raw counts, a
mean-variance trend, a leading axis tracking per-sample total signal (the
library-size warning) and [0,1]-bounded values, per view — so "what state is your
data in?" becomes something you read off the numbers rather than only ask about.

Per flag: say what was measured, name what it implies *for the model*, ask
whether that input state is intended. Then stop. Phrasing for each diagnostic is
in `references/data-prep.md`; read it before reporting.

What you are checking against — MOFA's assumptions about its input:

1. Per-sample scale effects aren't the dominant axis, or Factor 1 becomes a
   total-signal axis and everything subtler is crowded out behind it.
2. Each view's likelihood suits its scale. This one *is* yours to advise on —
   `scripts/recommend_likelihoods.R` measures what it turns on.
3. Technical structure has been handled however the user intends; MOFA has no
   batch term.
4. View widths within an order of magnitude, or `weight_views = TRUE` to
   compensate inside the model.
5. Features observed in no sample dropped; other missingness is fine.

If a check fails and the user asks you to fix it, implement what they choose
rather than a transform of your own.

## Stage 4 — Build the object

`create_mofa()` accepts a list of matrices (features × samples), a long
data.frame with `sample`/`feature`/`view`/`value` (+ optional `group`), a
`MultiAssayExperiment`, a `SingleCellExperiment`, or a `Seurat` object. Per-format
arguments and pitfalls: `references/input-formats.md`.

Then always `plot_data_overview()` and read it with the user. It puts sample
overlap and view dimensions in one picture, which is where dimensionality
imbalance and unexpectedly missing samples become obvious — while they are still
cheap to fix.

Name hygiene is enforced here and the errors are cryptic: view and group names
must not contain `/`; sample and feature names must be ASCII and unique within a
view.

For MEFISTO, attach covariates at this stage with `set_covariates()`.

## Stage 5 — Set options and train

Default to the **explicit four-step form**. It is what every published tutorial
uses, and it puts the option lists in variables the user can print, inspect and
modify — which is what makes a configuration reviewable:

```r
data_opts  <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)

model_opts$num_factors <- 15

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options     = data_opts,
                           model_options    = model_opts,
                           training_options = train_opts)
MOFAobject <- run_mofa(MOFAobject, outfile = "model.hdf5", use_basilisk = TRUE)
```

All five option groups with their real defaults, the `num_factors` reasoning, the
warnings `prepare_mofa()` raises on its own, and the compact `mofa2()`
alternative: `references/options.md`.

Two footguns worth carrying here, because they cost real time:

- **`use_basilisk` defaults to `FALSE`.** Training runs in Python (`mofapy2`)
  through reticulate, which is where most failures originate. Connection routes
  and failure triage: `references/troubleshooting.md`.
- **Always set `outfile`.** Training is the expensive step, and the default is a
  timestamped file in `tempdir()` that is easy to lose.

## Stage 6 — QC the model

Do this before interpreting anything, and don't skip it because training
succeeded — training always succeeds. Full diagnostic set and red-flag catalogue:
`references/model-qc.md`. The minimum:

```r
plot_variance_explained(model, x = "view", y = "factor")   # the central diagnostic
plot_variance_explained(model, plot_total = TRUE)          # total per view
plot_factor_cor(model)                                     # redundancy
```

Read the first plot with the user. A healthy model shows structure — some factors
shared across views, some private to one. Warning signs: a factor that dominates
a *single* view whose leading axis tracks per-sample totals (a depth effect in
that view); one view explaining nearly all variance (dimensionality imbalance); a
flat profile with no factor above a few percent (over-filtered data, too few
samples, or MOFA is genuinely the wrong tool here). A factor dominating *every*
view is not itself a warning sign — that is what a strong shared axis looks like,
biological or technical, so check what it correlates with before judging it.

## Stage 7 — Characterise

This is where you can be proactive. Attach metadata first, then loop per factor:

```r
samples_metadata(model) <- metadata_df   # must include a `sample` column
```

For each factor worth pursuing:

1. **Where does it live?** — variance explained across views and groups.
2. **What does it track?** — `correlate_factors_with_covariates()`,
   `plot_factor()` coloured by metadata.
3. **What drives it?** — `plot_top_weights()`, `plot_weights()` per view.
4. **Does the raw data agree?** — `plot_data_heatmap()`, `plot_data_scatter()`.
   This step is what separates a real factor from a plausible-looking artefact,
   so don't let it be skipped.
5. **What does it mean biologically?** — `run_enrichment()` →
   `references/enrichment.md`.

Function catalogue and interpretation guidance: `references/downstream.md`.

Two rules to state whenever you show factor or weight plots, because users
reliably over-read both. Factor *values* are meaningful only in relative terms —
the sign and the ordering of samples carry information, the absolute number does
not. And signs are arbitrary per run: a factor and its weights can flip wholesale
between runs without changing the solution, so signs may be compared within a
model but never across models.

## Reference files

Read these as the workflow reaches them, not upfront.

| File | When |
|---|---|
| `references/data-prep.md` | Stage 3 — what MOFA assumes of its input, how to report each diagnostic, likelihood choice |
| `references/input-formats.md` | Building the object from matrices / df / MAE / SCE / Seurat |
| `references/options.md` | All five option groups, real defaults, when to change them |
| `references/model-qc.md` | Diagnostics, red flags, robustness across seeds |
| `references/downstream.md` | Factor and weight characterisation, plots, export |
| `references/enrichment.md` | GSEA on weights; motif enrichment for ATAC |
| `references/mefisto.md` | Covariates, smoothness, sharedness, interpolation, warping |
| `references/troubleshooting.md` | reticulate/basilisk failures, GPU, scale-up |

Tutorials cited by name in these files are published, not local:
<https://biofam.github.io/MOFA2/tutorials.html>.

Scripts (run these to measure rather than infer; each has a `--help`):

| Script | When |
|---|---|
| `scripts/profile_mofa_data.R` | Stage 1 — dimensions, view balance, MEFISTO/GPU signals, a proposed configuration |
| `scripts/check_preprocessing.R` | Stage 3 — per-view readiness: raw counts, mean-variance trend, library-size axis, [0,1] values. Reports; does not prescribe |
| `scripts/recommend_likelihoods.R` | Stage 3 — per-view likelihood proposal from measured width/sparsity/range |