# A primer on what this project is doing

Written for team members who are not from computer science. It explains the data, the
methods and the words that keep coming up, in the order you need them. Nothing here
requires running code. If a term in a results file is unclear, it should be here.

The numbers for each dataset are in [dataset_facts.md](dataset_facts.md), regenerated
from the processed files by `scripts/dataset_facts.py`.

## 1. The data

A single-cell RNA sequencing experiment measures, for thousands of individual cells, how
much each gene is switched on. The result is a table:

- one **row per cell**, which is the sample or observation
- one **column per gene**, which is the feature or variable
- each entry is a **count**, roughly how many RNA molecules of that gene were seen in that cell

Two properties of this table drive everything else.

**It is wide.** Zeisel has 2,903 cells and 14,890 genes, so there are about five times more
variables than observations. Most classical statistics assumes the opposite.

**It is mostly zeros.** In PBMC14K about 94 percent of the entries are zero. A zero can mean
the gene is genuinely off in that cell, or that it was on but the experiment missed it. That
ambiguity is called **dropout**, and it matters for us because methods differ in how much they
rely on which genes are zero rather than on how large the non-zero values are.

**Counts per cell vary a lot.** One cell may have 2,000 total counts and another 60,000, mostly
for technical reasons. So we rescale before comparing cells.

## 2. What we do to the data before clustering

This is the preprocessing, in `src/scminer_clust/preprocess.py`. Four steps.

1. **Quality control.** Drop cells that look broken. We compute each cell's total count and
   number of detected genes, take logs, and drop a cell if either value is more than 4 median
   absolute deviations from the median. A median absolute deviation is a spread measure that is
   not thrown off by extreme values, which suits this data. On Zeisel this removes 18 of 3,005
   cells, all with very few detected genes.
2. **Normalization.** Rescale every cell to the same total, one million counts, which is called
   **counts per million** or CPM. Now a value means a share of that cell's output rather than
   a raw amount.
3. **Log transform.** Replace each value `x` with `log2(x + 1)`. Expression spans several orders
   of magnitude, and the log pulls that into a usable range. The `+ 1` keeps zeros at zero.
4. **Gene filtering.** Drop genes seen in fewer than 3 cells, since they carry no information
   for grouping cells.

**Highly variable genes (HVGs).** Many pipelines then keep only the 2,000 genes that vary most
across cells, on the view that the rest is noise. There are several ways to decide which genes
those are. We use Scanpy's three (`seurat`, `cell_ranger`, `seurat_v3`), and they disagree more
than you would expect. On Chung the `seurat` and `cell_ranger` sets share only 715 of 2,000
genes, and that changes the clustering result a lot. This is why we always say which gene set a
number came from.

## 3. Clustering, and the two methods we compare

**Clustering** means grouping the cells so that similar cells land together, without being told
the answer. The hope is that the groups correspond to **cell types**, such as B cells or
microglia.

Every clustering method needs a way to say how similar two cells are, and a way to turn those
similarities into groups.

### The standard pipeline: PCA plus Leiden

This is what most biologists run, and it is our baseline.

1. **PCA (principal component analysis)** compresses the thousands of genes down to about 50
   numbers per cell, keeping the directions along which cells differ most. It is a way to reduce
   noise and cost.
2. **Build a neighbor graph.** Join each cell to its 15 nearest cells in that compressed space.
   You now have a network of cells.
3. **Leiden** cuts that network into communities, which are groups with many links inside and
   few links outside.

Leiden takes a **resolution** parameter. Higher resolution gives more and smaller clusters. It
does not take the number of clusters directly, which matters below.

### scMINER's method: MICA

MICA is the clustering part of scMINER, the paper we are studying. It differs in how it measures
similarity between two cells.

Instead of PCA and straight-line distance, it uses **mutual information**. Take two cells. Sort
each cell's gene values into a small number of bins. Mutual information asks: if I tell you which
bin a gene falls into for cell A, how much does that tell you about its bin for cell B? If
knowing one tells you a lot about the other, the cells are similar. The appeal is that this
captures any consistent relationship, not just a straight-line one.

MICA then runs one of two pipelines depending on size:

- **Small-data mode**, up to 5,000 cells. Turn the distances into coordinates using **MDS
  (multidimensional scaling)**, which places cells in a space so the distances are preserved as
  well as possible, then run **k-means**. K-means needs to be told the number of clusters `k`.
- **Graph mode**, above 5,000 cells. Build a neighbor graph from the mutual information
  distances, learn coordinates from the graph with **node2vec**, then cut it with **Louvain**,
  which is Leiden's predecessor. This takes a resolution, like Leiden.

So MICA's knob is `k` on small data and a resolution on large data. This asymmetry is awkward
and we have to keep saying which mode a number came from.

## 4. The question we are actually asking

Here is the point of the project.

The scMINER paper reports that MICA beats the standard pipelines. But look at how its benchmark
was set up. Every method was either **given the true number of cell types**, or had its
resolution swept and only the runs that produced the true number were kept, with ties broken by
whichever matched the reference labels best.

That is a fair way to compare the clustering machinery with the number of clusters held fixed.
It is not what a biologist faces. Someone clustering new tissue does not know how many cell types
are there. That is usually the thing they are trying to find out.

**So our main question: when each method has to choose its own final clustering without the
reference labels, does MICA still beat PCA plus Leiden?**

The labels are still used, but only at the very end, to score the answer each method picked. They
are never used to make the choice. We report the paper's label-using setup separately, as an
upper bound on what each method could do.

A **label-free rule** is a recipe for making that choice without labels. We apply the same rule
to both methods. So far:

- **Default setting.** Just use whatever the software ships with. Leiden's resolution 1.0, MICA's
  1.822.
- **Best silhouette.** The **silhouette** of a cell compares how close it is to its own group
  against the nearest other group, from -1 to 1. Average it over cells, try every setting, keep
  the best. The catch is that it is computed in each method's own space, which may not be
  comparable.
- **Subsampling stability.** Rerun the method on many random 80 percent subsets. A setting is
  good if cells that grouped together keep grouping together. This follows a published tool
  called chooseR.
- **Transfer.** For each dataset, use the setting that worked best on the *other* datasets. No
  labels from the dataset being scored.

## 5. How we score a clustering

We compare the clusters a method produced against the reference labels.

**ARI (adjusted Rand index).** Take every pair of cells and ask whether the two groupings agree
about whether that pair belongs together. 1 means perfect agreement, 0 means no better than
chance. Chance correction is what the "adjusted" means, and it matters because splitting data
into many groups can look good by accident.

**AMI (adjusted mutual information).** The same idea built on information rather than pairs. We
report both because ARI reacts strongly to a method producing too many clusters, while AMI is
gentler about it. When the two disagree, that usually means the methods differ in how finely they
split, and that is worth saying out loud.

**A warning about the reference labels.** They are not always independent truth.

- **PBMC14K** labels came from physically sorting cells with antibodies before sequencing. That is
  as close to ground truth as this field gets.
- **Zeisel** labels came from the original authors clustering the same data. A method that
  resembles theirs will look better, which is a bias we have to mention.

## 6. Where each thing lives

| What | Where |
|---|---|
| The numbers for each dataset | `docs/dataset_facts.md` |
| One page per dataset, in our own words | `docs/dataset_cards/` |
| Every result ever produced, one row per run | `docs/results_tracker.csv` |
| Just the main comparison, for pivot tables | `docs/results_for_pivot.csv` |
| Tables ready for the paper | `latex/tables/` |
| Who does what, and by when | `WORKPLAN.md` |
| What has actually been done | `PROGRESS.md` |
| The current report | `latex/` |

## 7. Words you will meet in the result files

| Term | Meaning |
|---|---|
| ARI, AMI, NMI | How well a clustering matches the reference labels. Higher is better |
| silhouette | How well separated the clusters are, judged without labels |
| resolution | Leiden's and Louvain's knob. Higher gives more clusters |
| k | The number of clusters, when a method must be told it |
| HVG | Highly variable gene. The 2,000-gene subsets |
| `paper_*` | The input files the scMINER authors released |
| label-free | A choice made without using the reference labels |
| paper protocol | The paper's setup, where the labels pick the setting |
| graph mode, small-data mode | Which MICA pipeline was used |
| seed | The random starting point. Repeating with different seeds shows how much a result wobbles |
