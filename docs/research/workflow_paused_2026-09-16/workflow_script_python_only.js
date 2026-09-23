export const meta = {
  name: 'scminer-litreview-novelty',
  description: 'Deep literature review around scMINER, brainstorm novel publishable directions, adversarially verify novelty and feasibility',
  phases: [
    { title: 'Research', detail: '5 parallel literature/data sweeps' },
    { title: 'Ideate', detail: '3 lenses propose directions, then merge to a shortlist' },
    { title: 'Verify', detail: 'prior-art refuter + reviewer critic per idea' },
    { title: 'Synthesize', detail: 'rank, recommend, sketch plan and team split' },
  ],
}

const CONTEXT = `
PROJECT CONTEXT (read carefully):
- Graduate course "Clustering Algorithms", Missouri S&T, Fall 2026. Today is 2026-09-16. Semester ends mid December 2026.
- Team of 3. Utshab is the lead: CS background with strong information retrieval / RAG / LLM experience (ColBERT, TREC RAG, vLLM). Laraib and Rhea are from non-CS backgrounds. Only Utshab has compute: 2x NVIDIA RTX 6000 Ada (48 GB each), 120 CPU cores, 754 GB RAM, can run jobs about 12 hours uninterrupted. Keep things light: no training of big foundation models.
- Chosen paper: Pan Q, Ding L, Hladyshau S, et al. (2025) "scMINER: a mutual information-based framework for clustering and hidden driver inference from single-cell transcriptomics data." Nature Communications 16:4305. https://doi.org/10.1038/s41467-025-59620-6  (PDF: https://www.nature.com/articles/s41467-025-59620-6.pdf)
- Current proposal title: "Single-Cell Genomic Clustering Using scMINER with RAG-Assisted Biological Interpretation". The stated motivation is using cluster genes/markers to retrieve biomedical literature via RAG.
- Course rules: the paper must introduce/extend/adapt a clustering method. At least three datasets: Dataset 1 >= 5,000 samples and >= 20 variables (ordinary multivariate). Dataset 2 >= 2,000,000 samples and >= 10 variables (large-scale clustering). Dataset 3 >= 500 samples and >= 700 variables (high-dimensional clustering). Public or well documented datasets; document source, variables, preprocessing, labels.
- Goal: complete the course project AND attempt a publication if there is a genuinely novel, defensible contribution that is feasible by roughly Dec 2026 to early 2027.
RULES: This is READ-ONLY research. Do NOT create, edit or delete any files. Use WebSearch and WebFetch heavily (they are deferred tools: first call ToolSearch with query "select:WebSearch,WebFetch" to load them). Prefer primary sources (papers, bioRxiv/arXiv, GitHub READMEs/issues). Give a URL or DOI for every paper or claim you cite, with year and venue. Clearly mark anything you could not verify as UNVERIFIED. Cover 2023 through 2026 carefully, since the field moves fast.
`

const PY_ONLY = `
LANGUAGE CONSTRAINT (added after ideation): nobody on the team knows R. All team code must be Python. scMINER's clustering (MICA) and network tool (SJARACNe) are already Python/C++. Its R-only steps (QC, CPM+log2 normalization, driver activity as mean of targets, t-test for hidden drivers) must be reimplemented in Python. R-only tools (Seurat, CHOIR, sc-SHC, recall, chooseR, MultiK, NetBID2) may appear only as optional black-box baselines run through a Python-launched Rscript wrapper that Utshab does not have to maintain, or be replaced by Python equivalents (Scanpy, rapids-singlecell, SCCAF). Penalize directions whose core method or essential baselines require writing or maintaining R code.
`

const RESEARCH = [
  { key: 'scminer-deep-read', prompt: `${CONTEXT}
TASK: Deep read of scMINER and its software. Fetch the paper (and supplementary info if reachable) plus GitHub repos (jyyulab/scMINER, jyyulab/MICA, SJARACNe, NetBID2, and the scMINER docs site). Report precisely:
1. MICA clustering pipeline step by step: preprocessing, how mutual information between cells is estimated (binning scheme, number of bins, on which genes), MI to distance transform, dimension reduction paths (MDS for small data vs graph embedding / node2vec on kNN graph for large data, and the cell-count threshold), clustering algorithm (k-means / Louvain / Leiden), how k or resolution is chosen, any consensus steps, defaults.
2. SJARACNe network inference and NetBID activity / hidden driver inference: inputs, outputs, how per-cluster networks are built, compute cost.
3. Computational complexity, reported runtime and memory, the largest dataset they ran, any downsampling used for large data, GPU support or not. Would it run on 2,000,000 cells? Why or why not? Be concrete.
4. Benchmark datasets used (names, cell counts, gene counts, label sources, accessions), baselines compared (Seurat, Scanpy, SC3, etc.), metrics (ARI, NMI, etc.), and headline numbers.
5. Software stack: languages, install method (pip/conda/R), Python and R versions, known install issues from GitHub issues, input formats (h5ad, SparseEset).
6. Limitations stated in the Discussion and any future work the authors mention. Limitations you infer (for example handling of dropout zeros in binned MI, quadratic cost, resolution choice).
7. Papers from 2025-2026 that cite or extend scMINER or MICA (use Google Scholar / Semantic Scholar / Europe PMC "cited by" searches).
Return a detailed, well organized report (1500-3000 words) with citations.` },

  { key: 'scalable-mi-clustering', prompt: `${CONTEXT}
TASK: Literature sweep on (a) information-theoretic and mutual-information based similarity/clustering and (b) scalable clustering for million-cell single-cell data.
Cover at least:
- MI or other information-theoretic cell-cell or gene-gene similarity for scRNA-seq (MICA, scMINER, proportionality, "impact of similarity metrics on scRNA-seq clustering" benchmarks such as Kim et al. 2019 and Skinnider et al. 2019, and any 2022-2026 benchmarks of distance metrics). Also MI estimators for sparse count data (binning, KSG kNN estimator, zero-inflation aware MI, adaptive partitioning), and fast/approximate/GPU MI computation.
- Scalable clustering at atlas scale: sketching and subsampling (geosketch, Hopper, Kernel Herding sketching, scSampler), landmark / Nystrom / anchor-graph spectral clustering, approximate kNN (FAISS, HNSW, pynndescent), GPU single-cell stacks (rapids-singlecell, cuML), Leiden on huge graphs, mini-batch k-means, SCALEX-like approaches, and how atlas projects cluster 1M-100M cells (CELLxGENE Census, Tabula Sapiens, Allen Brain Cell atlas).
- Any prior work that already scales mutual-information distances to millions of samples (in any field, not just biology), for example random projection MI, locality sensitive hashing for MI, sketching for MI matrices, landmark MI.
- Cluster-then-propagate strategies: cluster a sketch with an expensive metric, then assign remaining points; theoretical or empirical guarantees on how rare populations survive sketching.
End with a GAP ANALYSIS: what specifically is NOT yet done, with evidence of the search you did (queries tried). Be precise and skeptical. 1500-2500 words with citations.` },

  { key: 'llm-rag-annotation', prompt: `${CONTEXT}
TASK: Literature sweep on LLM and retrieval-augmented methods for interpreting clusters / annotating cell types, 2023-2026.
Cover at least: GPTCelltype (Hou & Ji, Nature Methods 2024), GeneAgent (Nature Methods 2025), CellAgent, CASSIA, scExtract, CellTypeAgent, mLLMCelltype, LICT, AnnDictionary, scAgent, CellVoyager, Biomni, BioDiscoveryAgent, Cell2Sentence / C2S-Scale, scGPT / Geneformer annotation, GenePT / gene embedding approaches, retrieval over CellMarker 2.0 / PanglaoDB / PubMed, and any 2025-2026 benchmarks of LLM cell type annotation (including hallucination and grounding evaluation).
Answer these specific questions with evidence:
Q1. Does any method use regulatory network activity or inferred driver / transcription factor activity (e.g. from SCENIC, VIPER, DoRothEA, NetBID) as the query/evidence for LLM or RAG annotation, instead of or in addition to differentially expressed marker genes?
Q2. Does any method feed the LLM/RAG annotation verdict back into the clustering itself (choose resolution, merge or split clusters, detect over-clustering)?
Q3. Does any method use retrieval to verify or cite its annotations, and how is faithfulness measured?
Q4. What datasets, label ontologies (Cell Ontology), and metrics do they use? What are typical accuracy numbers?
Q5. Which open-weight LLMs (Qwen, Llama, etc.) have been evaluated, and how do they compare to GPT-4 class models? Could a local 30B-70B model on 2x48 GB GPUs be competitive?
End with a GAP ANALYSIS listing what is not yet done, with the queries you tried. 1500-2500 words with citations.` },

  { key: 'datasets', prompt: `${CONTEXT}
TASK: Find concrete, downloadable public datasets for the three course regimes, favoring single-cell transcriptomics so the project stays coherent, and favoring datasets with trustworthy cell type labels.
For each candidate give: name, paper and year, number of samples (cells) and variables (genes), label type and granularity, source/accession, direct download format (h5ad preferred, e.g. CELLxGENE, figshare, GEO, 10x website, scvi-tools datasets, Zenodo), file size, license, whether scMINER or common clustering benchmarks used it, and preprocessing pitfalls.
- Dataset 1 (>= 5,000 cells, >= 20 variables, ordinary): e.g. 10x PBMC sets, Zheng 2017 sorted PBMC, Tabula Muris subsets, human pancreas.
- Dataset 2 (>= 2,000,000 cells, >= 10 variables, large-scale): e.g. Cao 2019 mouse organogenesis (MOCA, ~2.06M cells), Cao 2020 human fetal atlas (~4M), CELLxGENE Census slices, Tabula Sapiens v2, Allen Brain Cell atlas (ABC), Parse Biosciences 10M PBMC, Tahoe-100M, HLCA core/extended. Check which have high quality labels and which can be downloaded as one file.
- Dataset 3 (>= 500 cells, >= 700 variables, high-dimensional): classic gold-standard sets like Usoskin, Zeisel, Kolodziejczyk, Baron, Segerstolpe, Muraro, Darmanis, Chu, plus scMINER benchmark sets.
Also: estimate memory to hold 2M cells x 2,000 HVGs as float32 dense and as sparse, and whether it fits in 754 GB RAM and in 48 GB GPU memory.
Also say whether a non-biology dataset would be wiser for Dataset 1 or Dataset 2 given the course wording, and list 2 strong non-biology alternatives for each regime just in case (with sizes and labels).
Recommend a primary and a backup pick per regime with reasons. 1200-2500 words with URLs.` },

  { key: 'resolution-validity', prompt: `${CONTEXT}
TASK: Literature sweep on choosing clustering resolution / number of clusters and judging cluster validity in single-cell data, plus evaluation pitfalls.
Cover at least: CHOIR (Nature Genetics 2025), scSHC significance of hierarchical clustering (Grabski et al. 2023), chooseR, clustree, MultiK, scCCESS, SC3 consensus, IKAP, "clustering deserves more attention" style critiques, double-dipping / post-selection inference for clusters (count splitting, Neufeld et al.), cluster stability via subsampling, silhouette and its problems in high dimension, ARI/NMI/AMI pitfalls with imbalanced clusters, rare-cell-type detection metrics, and label hierarchy aware evaluation (Cell Ontology distance).
Also look for knowledge-guided or annotation-guided resolution selection: methods that use marker gene databases, gene set enrichment, or LLMs to decide whether clusters should be merged or split.
And look for benchmark papers of clustering methods in 2023-2026 (e.g. large scRNA-seq clustering benchmarks) and where MI-based methods rank if included.
End with a GAP ANALYSIS with the queries you tried. 1200-2500 words with citations.` },
]

phase('Research')
const research = await parallel(RESEARCH.map(r => () =>
  agent(r.prompt, { label: `research:${r.key}`, phase: 'Research' }).then(t => ({ key: r.key, text: t }))
))
const researchOk = research.filter(Boolean).filter(r => r.text)
log(`Research done: ${researchOk.length}/${RESEARCH.length} reports`)
const RESEARCH_DUMP = researchOk.map(r => `\n\n===== REPORT: ${r.key} =====\n${r.text}`).join('')

const IDEA_SCHEMA = {
  type: 'object',
  properties: {
    ideas: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          title: { type: 'string' },
          one_line: { type: 'string' },
          problem_and_gap: { type: 'string', description: 'What real problem, what evidence it is unsolved' },
          method_sketch: { type: 'string', description: 'Concrete method, how it extends or adapts scMINER/MICA' },
          experiments: { type: 'string', description: 'Datasets per regime, baselines, metrics, ablations' },
          compute_estimate: { type: 'string' },
          team_fit: { type: 'string', description: 'Which parts non-CS teammates can own' },
          novelty_claim: { type: 'string', description: 'The precise claim a reviewer would check' },
          closest_prior_work: { type: 'string' },
          risks: { type: 'string' },
          venue_ideas: { type: 'string' },
        },
        required: ['title', 'one_line', 'problem_and_gap', 'method_sketch', 'experiments', 'compute_estimate', 'team_fit', 'novelty_claim', 'closest_prior_work', 'risks', 'venue_ideas'],
      },
    },
  },
  required: ['ideas'],
}

const LENSES = [
  { key: 'clustering-method-first', text: 'You are a senior clustering / machine learning researcher. Prioritize a genuine, defensible contribution to the CLUSTERING METHOD itself (e.g. making MI based clustering scale to millions of cells, better MI estimation for sparse counts, principled resolution selection, stability). The course grades clustering, so the core novelty should be a clustering contribution.' },
  { key: 'interpretation-rag-first', text: 'You are a senior researcher at the intersection of information retrieval, LLMs, and computational biology. Utshab has strong RAG / IR experience. Prioritize ideas where retrieval or LLM grounding creates a genuinely new capability tied to clustering (not yet another LLM cell type annotator), while keeping a clear clustering contribution.' },
  { key: 'feasibility-first', text: 'You are a pragmatic PI who has shipped many student papers. Prioritize ideas that a 3 person team with 2 non-CS members can finish by mid December 2026 on one workstation, that satisfy all three dataset regimes naturally, and that still have a real shot at a workshop or journal publication. Avoid anything needing heavy model training.' },
]

phase('Ideate')
const ideaSets = await parallel(LENSES.map(l => () =>
  agent(`${CONTEXT}
${l.text}

Below are five research reports produced by colleagues. Read all of them closely. Then act like an expert problem finder: look for real, specific pain points and unexplored gaps (not generic "apply X to Y"). Propose 4 distinct research directions that build on scMINER. Each must satisfy the course rules (extend/adapt a clustering method, three dataset regimes) and be light on compute. You may do a few extra WebSearch checks to sharpen ideas, but do not write files. Be concrete and honest about novelty risk.
${RESEARCH_DUMP}`, { label: `ideate:${l.key}`, phase: 'Ideate', schema: IDEA_SCHEMA })
    .then(r => r ? r.ideas.map(i => ({ ...i, lens: l.key })) : [])
))
const allIdeas = ideaSets.filter(Boolean).flat()
log(`Ideation produced ${allIdeas.length} ideas`)

const SHORTLIST_SCHEMA = {
  type: 'object',
  properties: {
    candidates: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          id: { type: 'string' },
          title: { type: 'string' },
          merged_from: { type: 'string' },
          full_description: { type: 'string', description: 'Self-contained: problem, gap, method, experiments per regime, compute, team fit, novelty claim, closest prior work, risks' },
          novelty_claim: { type: 'string' },
        },
        required: ['id', 'title', 'merged_from', 'full_description', 'novelty_claim'],
      },
    },
    dropped: { type: 'string', description: 'Ideas dropped and why' },
  },
  required: ['candidates', 'dropped'],
}

const shortlist = await agent(`${CONTEXT}
${PY_ONLY}
You are the program committee chair. Below are ${allIdeas.length} proposed research directions from three different lenses. Merge duplicates and near duplicates, combine complementary ideas where that makes a stronger single project, and select the 5 strongest candidates for independent verification. Balance: at least one must have its core novelty in the clustering method itself, at least one should use Utshab's RAG strength, and all must be feasible for this team. Each candidate description must be fully self-contained because verifiers will not see anything else. List what you dropped and why.

IDEAS JSON:
${JSON.stringify(allIdeas, null, 1)}`, { label: 'ideate:shortlist', phase: 'Ideate', schema: SHORTLIST_SCHEMA })

const candidates = (shortlist && shortlist.candidates) ? shortlist.candidates.slice(0, 5) : []
if (shortlist && shortlist.candidates && shortlist.candidates.length > 5) log(`Shortlist had ${shortlist.candidates.length}, verifying top 5 only`)

const PRIOR_ART_SCHEMA = {
  type: 'object',
  properties: {
    refuted: { type: 'boolean', description: 'true if prior work already does essentially this' },
    novelty_score: { type: 'integer', description: '1 (already done) to 5 (clearly new)' },
    closest_works: { type: 'array', items: { type: 'object', properties: { citation: { type: 'string' }, url: { type: 'string' }, overlap: { type: 'string' } }, required: ['citation', 'url', 'overlap'] } },
    queries_tried: { type: 'array', items: { type: 'string' } },
    what_remains_novel: { type: 'string' },
    how_to_sharpen_claim: { type: 'string' },
  },
  required: ['refuted', 'novelty_score', 'closest_works', 'queries_tried', 'what_remains_novel', 'how_to_sharpen_claim'],
}
const CRITIC_SCHEMA = {
  type: 'object',
  properties: {
    feasibility_score: { type: 'integer', description: '1-5 for this team, compute and a mid December 2026 course deadline' },
    significance_score: { type: 'integer', description: '1-5, would reviewers at a good venue care' },
    course_fit_score: { type: 'integer', description: '1-5, fits clustering course rules and three regimes naturally' },
    fatal_flaws: { type: 'string' },
    compute_reality_check: { type: 'string', description: 'Concrete estimate of runtimes and memory on 2x48GB GPU, 120 cores, 754GB RAM; which step is the bottleneck' },
    tasks_for_non_cs_members: { type: 'string' },
    minimum_viable_paper: { type: 'string', description: 'Smallest version that is still publishable' },
    realistic_venues: { type: 'string', description: 'Named venues with approximate 2026-2027 deadlines, marked UNVERIFIED if unsure' },
    reviewer_objections: { type: 'string' },
  },
  required: ['feasibility_score', 'significance_score', 'course_fit_score', 'fatal_flaws', 'compute_reality_check', 'tasks_for_non_cs_members', 'minimum_viable_paper', 'realistic_venues', 'reviewer_objections'],
}

phase('Verify')
const verified = await pipeline(candidates,
  c => parallel([
    () => agent(`${CONTEXT}
You are an adversarial prior-art hunter. Your job is to REFUTE the novelty of the research direction below. Search hard (at least 12 distinct queries across Google Scholar style searches, bioRxiv, arXiv, PubMed/Europe PMC, Semantic Scholar, GitHub). Look for work from 2020-2026 in single-cell AND in other fields (general ML, IR, data mining) that already does the core idea. Default to refuted=true only if you find work that does essentially the same thing. Otherwise state exactly what remains novel and how to phrase the claim so it survives review.

CANDIDATE:
Title: ${c.title}
Novelty claim: ${c.novelty_claim}
${c.full_description}`, { label: `prior-art:${c.id}`, phase: 'Verify', schema: PRIOR_ART_SCHEMA }),
    () => agent(`${CONTEXT}
${PY_ONLY}
You are a tough but fair reviewer and a pragmatic advisor. Critique the research direction below for feasibility with this exact team and hardware, significance, and course fit. Do concrete compute math (cell counts, genes, pairwise costs, GPU memory). Name fatal flaws plainly. You may use WebSearch to check tool capabilities, dataset sizes, and venue deadlines.

CANDIDATE:
Title: ${c.title}
Novelty claim: ${c.novelty_claim}
${c.full_description}`, { label: `critic:${c.id}`, phase: 'Verify', schema: CRITIC_SCHEMA }),
  ]).then(([pa, cr]) => ({ candidate: c, prior_art: pa, critic: cr }))
)
const verifiedOk = verified.filter(Boolean)
log(`Verified ${verifiedOk.length} candidates`)

phase('Synthesize')
const synthesis = await agent(`${CONTEXT}
${PY_ONLY}
You are the senior advisor writing the final research brief for Utshab (the team lead). Use the research reports, the shortlist, and the verification results below. Produce a detailed markdown brief with these sections:
1. Executive summary (what scMINER does, where it falls short, the recommended direction in 3-4 sentences).
2. Literature review: organized by theme (MI based clustering and scMINER internals, scalable clustering, LLM/RAG annotation, resolution selection and validity). Keep real citations with URLs/DOIs. Mark UNVERIFIED items.
3. Ranked table of the verified candidates with novelty, feasibility, significance, course fit scores, and a one line verdict each.
4. RECOMMENDED DIRECTION in depth: final framing, precise novelty claim that survives the prior-art check, method, experiments per dataset regime (named datasets with sizes), baselines, metrics, ablations, compute plan with runtime estimates, risks and fallbacks, minimum viable paper, target venues. Also a BACKUP direction in brief.
5. Implementation plan skeleton: phases/milestones from mid September to mid December 2026, with what code runs where.
6. Team work split: concrete tasks for Utshab (all GPU/compute runs, core code), and for two non-CS teammates (Laraib, Rhea). Non-CS tasks should be real and valuable (dataset documentation, biological validation of clusters against marker databases, literature tables, writing sections, figure review, manual annotation audits), with the skills they need and how to hand off work without them running GPU code.
7. Open questions the lead must decide.
Be honest where novelty is thin. Do not pad.

RESEARCH REPORTS:${RESEARCH_DUMP}

SHORTLIST (dropped notes: ${shortlist ? shortlist.dropped : 'n/a'}):
${JSON.stringify(verifiedOk, null, 1)}`, { label: 'synthesize:brief', phase: 'Synthesize' })

return {
  synthesis,
  verified: verifiedOk.map(v => ({
    id: v.candidate.id, title: v.candidate.title,
    refuted: v.prior_art ? v.prior_art.refuted : null,
    novelty: v.prior_art ? v.prior_art.novelty_score : null,
    feasibility: v.critic ? v.critic.feasibility_score : null,
    significance: v.critic ? v.critic.significance_score : null,
    course_fit: v.critic ? v.critic.course_fit_score : null,
  })),
  research_keys: researchOk.map(r => r.key),
}
