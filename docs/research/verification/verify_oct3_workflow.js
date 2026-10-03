export const meta = {
  name: 'verify-oct3-claims',
  description: 'Recompute, attack and code-review the Oct 3 results before they are written up',
  phases: [{ title: 'Verify', detail: 'numbers, skeptic, code review' }],
}
const REPO = '/home/ugdf8/coursework/26Fall/Clustering_Algorithms/scminer-clustering'
const CLAIMS = '/tmp/claude-1599031162/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/scratchpad/claims_oct3.md'
const CTX = `Read ${CLAIMS} first. Repo ${REPO}. Python: cd ${REPO} && source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer && OMP_NUM_THREADS=4 python ...  (MICA's own env is mica-legacy).
READ-ONLY: do not modify, create or delete files in the repo or data folders, and do not start or stop processes. Other jobs are running on this shared machine: keep your own computations light (cap threads at 4, avoid recomputing distance matrices).`
const SCHEMA = { type: 'object', properties: {
  claims: { type: 'array', items: { type: 'object', properties: {
    id: { type: 'string' }, verdict: { type: 'string', enum: ['supported', 'supported with correction', 'not supported', 'unclear'] },
    evidence: { type: 'string' }, correction: { type: 'string' } }, required: ['id', 'verdict', 'evidence', 'correction'] } },
  other_problems: { type: 'array', items: { type: 'string' } }, summary: { type: 'string' } },
  required: ['claims', 'other_problems', 'summary'] }
phase('Verify')
const out = await parallel([
  () => agent(`${CTX}
You are an independent numbers checker. Recompute every number in claims R1 to R4, B1 to B4, C1 to C3 and D1 to D4 from the raw files (results/*/metrics.csv, results/uncertainty/*.csv, MICA output folders), not from the summaries where avoidable. Show the arithmetic for the means, counts and ranges. Flag any stale rows (old experiment folders that no longer match the current code, for example MICA rows scored before the 100-restart change, or label-free rows from before n_init changes).`, { label: 'verify:numbers', phase: 'Verify', schema: SCHEMA }),
  () => agent(`${CTX}
You are an adversarial skeptic. Try to refute or weaken each claim with alternative explanations. In particular: is the 100-restart choice for MICA fair when the k-means baseline uses 10 restarts and Leiden picks the best of 200 resolutions (does it favor MICA, and do the paper's published numbers really correspond to the best fit or could that be coincidence)? Is Leiden's failure on the small datasets a property of Leiden or of our settings (n_neighbors=15, 50 principal components, scaling, on 124 to 317 cells)? Is silhouette a fair label-free rule across embeddings of different kinds (MDS, node2vec, PCA), and does the k=2..20 range bias it? Do label sources (Zeisel, Usoskin, Buettner, Pollen and others from clustering, versus design or sorting) explain any pattern? Are the bootstrap intervals and Wilcoxon tests used correctly (what do they and don't they capture)? For each claim give the strongest honest wording.`, { label: 'verify:skeptic', phase: 'Verify', schema: SCHEMA }),
  () => agent(`${CTX}
You are a code reviewer. Review the code written on Oct 2 and 3 for bugs that could make these numbers wrong: src/scminer_clust/mi.py, src/scminer_clust/diagnostics.py, src/scminer_clust/selection.py, scripts/run_label_free.py, scripts/mica_label_free.py, scripts/mica_seeded_kmeans.py, scripts/bootstrap_gaps.py, scripts/run_diagnostics.py, scripts/summarize_label_free.py, scripts/summarize_gene_sets.py, data.load_paper_input in src/scminer_clust/data.py, scripts/prepare_paper_inputs.py, scripts/run_baselines.py. Check especially: cell alignment between partitions and labels everywhere (bootstrap, label-free, node2vec embedding row order), that mica_label_free sorts the node2vec embedding into the same cell order as the clustering files, that bootstrap pairs the right MICA and Leiden partitions for each gene setting, that the detection-only distance and the co-detection share are computed as described, that thinning uses raw counts and the same genes, and the log base and counts layer for the authors' inputs. Use small in-memory checks to confirm suspected bugs.`, { label: 'verify:code', phase: 'Verify', schema: SCHEMA }),
])
return out.filter(Boolean)
