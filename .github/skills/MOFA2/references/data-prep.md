# Data readiness

Preprocessing decides more about a MOFA result than any model option does. The
model cannot tell a technical axis from a biological one — it decomposes whatever
variance it is given, in proportion to how much there is.

**Scope.** Check the input and report what the numbers show; don't choose how the
data should be transformed. Which normalisation, variance-stabilisation,
batch-correction or feature-selection method fits depends on the assay, the
platform, the design and processing that already happened upstream — the
analyst's knowledge, not something recoverable from a matrix. So run the
diagnostics, state what they found, name the MOFA-specific consequence ("an
unremoved depth effect will show up as Factor 1"), and ask what state each view
is in. Implementing a step the user has chosen is fine; proposing one unprompted
is not.

Two scripts make this measurement rather than inference from a description:

- **`scripts/check_preprocessing.R`** — flags raw counts, a mean-variance trend,
  a leading axis tracking per-sample total signal (the library-size warning sign)
  and [0,1]-bounded values. Run it on the input before building anything.
- **`scripts/recommend_likelihoods.R`** — measures each view's width, sparsity
  and range, which is what the likelihood decision turns on.

Both print the measurements behind every call. Use them to ground the
conversation; the analyst owns the decisions.

## What MOFA assumes about its input

Properties of the model, not preferences. Check each, report the ones that fail,
leave the remedy to the user.

1. **Per-sample scale effects are not the dominant axis.** If library size or
   depth still drives the data, Factor 1 becomes a total-signal axis and subtler
   variation is crowded out behind it — the most common cause of a meaningless
   Factor 1. (MOFA2 FAQ, "How do I normalise the data?")
2. **Each view's likelihood matches its scale** — see below. For sequencing
   counts, preprocessing exists to make the gaussian likelihood and its
   homoscedasticity assumption appropriate (MEFISTO preprocessing vignette). For
   bounded values, report the bound and leave the judgement to the analyst.
3. **Technical structure has been handled however the analyst intends.** MOFA has
   no batch term, so anything left in the data can become a factor. If they do
   plan to remove it, one ordering fact is MOFA-relevant and easy to get
   backwards: removal has to precede feature selection, since an unremoved batch
   effect distorts which features look variable and bakes batch structure into
   the feature set. (MOFA2 FAQ.)
4. **View widths within an order of magnitude.** Wider modalities are
   over-represented in the factors: a 20,000-feature RNA view beside a
   50-feature mutation view yields factors that are essentially about RNA.
   (MOFA2 FAQ, "My data sets have different dimensionalities.")
5. **Missing values are fine, empty features are not** — see below.

## Reporting what the diagnostics found

Translate each signal into its consequence for the model, and stop there.

| Signal | What it means for MOFA | What to say |
|---|---|---|
| Non-negative integers, wide range | Looks like raw counts; gaussian fits them poorly and depth effects are still present | Report it, ask whether this is the intended input state |
| Leading axis tracks per-sample total signal | Factor 1 will likely be that axis | Flag it — but note it may be real biology where the per-sample total is itself meaningful (global drug sensitivity, total abundance) |
| Strong mean-variance trend | Conflicts with gaussian's homoscedasticity assumption | Report the correlation, let the user judge whether it matters |
| Values bounded in [0,1] | Gaussian has unbounded support | Report the bound; the user knows whether these are rates, proportions, or already appropriate |
| Continuous with negative values | Consistent with an already-transformed view | Note it, confirm what the transform was |
| Binary | Not a scale question but a likelihood one | → next section |

Where the numbers leave a view ambiguous, ask. "Counts", "log-CPM", "VST", "beta
values", "CLR-transformed" behave differently in the model, and "RNA-seq data"
alone does not distinguish them.

## Choosing a likelihood

This one is yours to advise on — it is a choice about how the model reads a view,
not about how the view should be transformed.

MOFA and MEFISTO support **gaussian** (continuous), **poisson** (counts) and
**bernoulli** (binary). The non-gaussian ones use a Gaussian approximation to keep
variational inference fast in the non-conjugate setting (Seeger & Bouchard; noted
in the MEFISTO preprocessing vignette) — a mechanism, not a defect.

Per view:

- **Continuous → gaussian.** The default, and where most views end up.
- **Narrow, biologically binary (somatic mutations) → bernoulli.** The CLL
  tutorial does exactly this: `model_opts$likelihoods["Mutations"] <-
  "bernoulli"`. Sparsity is not an argument against it — rare 0/1 events are what
  bernoulli is for.
- **Very wide binary (scATAC peaks)** — bernoulli is valid but slow across tens
  of thousands of features. Report that cost. Gaussian applies to a continuous
  version of the view; the 10x Multiome tutorial TF-IDF-normalises peaks for
  exactly this reason, but whether that suits the data is the user's call.
- **Counts** — poisson exists, and the tutorials rarely use it on real data (only
  a simulated MEFISTO demo); their usual route for sequencing counts is library
  size correction plus variance stabilisation, then gaussian. Report what the
  view looks like and let the user choose which input to supply.

`recommend_likelihoods.R` measures the width, sparsity and range these turn on, so
narrow binary and wide sparse views — indistinguishable to a naive check — are
separated by numbers rather than by the modality's name.

> **Note on the package default.** Automatic likelihood inference is currently
> **disabled** — the `.infer_likelihoods()` call is commented out in
> `prepare_mofa.R`, so every view defaults to gaussian regardless of its content.
> Any non-gaussian view must be set explicitly in `model_opts$likelihoods`. Some
> tutorials still say likelihoods are "inferred automatically"; that is not true
> of the current code.

## View balance

Width imbalance is a MOFA property rather than a preprocessing preference, so
surface it explicitly. Two levers exist, both the user's to pull:

- Reduce the wide view (feature selection, or splitting it into narrower,
  biologically distinct views). Which features to keep is data-specific — report
  the widths and the ratio, not a target number.
- `weight_views = TRUE` in the data options, which rebalances inside the model
  rather than in the data. → `references/options.md`.

One caveat specific to **multi-group** models: the framework centres features
within each group before fitting, so feature selection on uncentred data tends to
pick features for between-group differences the model then removes — features
chosen for variance that no longer exists. (MOFA2 FAQ.) State this whenever
groups and feature selection are both in play.

## Missing values

MOFA drops missing values from the likelihood — no hidden imputation step, and
matrix factorisation is robust to substantial missingness (MOFA2 FAQ, "Does MOFA
handle missing values?"). Two things still need attention:

- Features with no observation in *any* sample contribute nothing, trigger
  warnings, and should be dropped.
- Samples missing an entire view are fine, as long as a good fraction of samples
  have all views.

How the missingness arose is worth asking about: MOFA treats `NA` as
uninformative, which suits missing-at-random better than values missing because
they fell below a detection limit. Raise the distinction; the user decides what
to do with it.