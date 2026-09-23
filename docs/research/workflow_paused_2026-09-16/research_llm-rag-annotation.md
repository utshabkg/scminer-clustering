**LLM and retrieval-augmented cluster interpretation and cell type annotation, 2023 to 2026 (sweep done 2026-09-16)**

**How I checked the sources.** Nature, Springer and Genome Biology blocked automated fetches, so several numbers come from the PMC or bioRxiv copies of the same papers. Numbers in a published version may differ. The automatic page summaries were unreliable in at least one case. For CellMaster the summary said the tool uses RAG, and the PDF does not support that. So I read the full PDFs for CellMaster, GATHER and HypoGeneAgent myself. The web search budget ran out near the end, so a few items stay UNVERIFIED and are listed at the bottom.

## 1. Main methods found

**Prompting one LLM with a marker list**
- **GPTCelltype.** Hou and Ji, *Nature Methods* 21:1462 (2024), doi:10.1038/s41592-024-02235-4 (full text: PMC11310073).
  - GPT-4 gets the top 10 Wilcoxon DEGs for each cluster. Tested on 10 datasets across 5 species.
  - Labels fully or partially matched the manual labels for over 75% of cell types in most studies. Scoring was full match 1, partial 0.5, mismatch 0, using Cell Ontology (CL) names.
  - It gave the same answer for the same input 85% of the time. The paper warns about hallucination and about the unknown training corpus.
- **DeepCellSeek benchmark.** Xiao et al., *Briefings in Bioinformatics* 26(6):bbaf677 (2025), https://academic.oup.com/bib/article/26/6/bbaf677/8380197.
  - 11 API models on 34 datasets. Perfect match on broad types: Kimi-K2 83.5%, GPT-5 80.6%. SingleR, ScType and CellMarker2.0 were all under 20%.
  - The top 30 genes ranked by p-value plus log2FC beat the top 10. Broad types held up with 60% marker contamination. Subtypes did not.

**Consensus across several LLMs**
- **LICT.** *Communications Biology* (2025), doi:10.1038/s42003-025-08745-x (preprint arXiv:2409.15678).
  - Uses GPT-4, Claude 3 Opus, Gemini 1.5 Pro, ERNIE 4.0 and Llama-3-70B, plus a "talk-to-machine" retry loop.
  - Credibility rule: a label counts as reliable if more than four of the markers the LLM proposed are expressed in at least 80% of the cluster's cells.
- **mLLMCelltype.** bioRxiv (2025), doi:10.1101/2025.04.10.647852.
  - The models deliberate over several rounds. It scored 77.3% against 61.3% for GPTCelltype across 50 datasets.
  - Clusters with low agreement or high entropy are flagged. It has no RAG. The authors list RAG as future work.
- **AnnDictionary.** *Nature Communications* 16:9511 (2025), doi:10.1038/s41467-025-64511-x (preprint bioRxiv 10.1101/2024.10.10.617605).
  - 15 models on Tabula Sapiens v2, with Leiden resolution fixed at 0.5.
  - Claude 3.5 Sonnet reached 84.0% agreement (kappa 0.72). Agreement was judged by an LLM that was never checked against humans. Accuracy tracked model size.

**Multi-agent pipelines**
- **CASSIA.** *Nature Communications* (2025), doi:10.1038/s41467-025-67084-x (preprint bioRxiv 10.1101/2024.12.04.626476).
  - Agents: annotator, validator, scorer and reporter. Optional agents: Annotation Boost, Subclustering, and a RAG agent that queries CellMarker2.0 and ontologies (according to the preprint).
  - Scores are hierarchical on CL, covering 970+ cell types.
- **CellTypeAgent.** arXiv:2505.08844 (2025). The LLM proposes candidate types. The final label is the candidate whose markers have the highest mean expression in CZ CELLxGENE. Tested on 9 datasets, 303 cell types and 36 tissues.
- **CyteType.** bioRxiv (2025), doi:10.1101/2025.11.06.686964. Competing hypotheses are checked against CL, PubMed and marker databases, with pathway enrichment. It reports confidence scores and heterogeneity flags. Details are only partly verified.
- **scExtract.** *Genome Biology* (2025), doi:10.1186/s13059-025-03639-x. It reads the source paper to guide processing and annotation, then uses the labels as priors for integration (scanorama-prior, cellhint-prior).
- **GPTBioInsightor.** *Bioinformatics Advances* 6(1):vbag025 (2026), PMC12975716. The LLM input is DEGs plus pathway enrichment computed with GSEApy.
- **scPilot.** arXiv:2602.11609 (2026). Annotation, trajectory and TF-target tasks on its own scBench.
- **CellMaster.** arXiv:2602.13346 (2026). Detailed under Q2.
- **General agents.**
  - CellAgent: arXiv:2407.09811, ICLR 2026 poster.
  - Biomni: bioRxiv 10.1101/2025.05.30.656746 (2025). Its eight test scenarios include single-cell annotation.
  - CellVoyager: *Nature Methods* 23:749 (2026), doi:10.1038/s41592-026-03029-6. It does open-ended reanalysis, not annotation.
  - BioDiscoveryAgent: ICLR 2025, arXiv:2405.17631. It designs perturbation experiments with PubMed search and does not annotate cells.

**Methods grounded in retrieval or knowledge graphs**
- **GeneAgent.** *Nature Methods* 22 (2025), doi:10.1038/s41592-025-02748-6, PMC12328209. It explains gene sets, not cell types. It checks its own claims against 18 databases.
- **scRAG.** *Findings of ACL* 2025, https://aclanthology.org/2025.findings-acl.53/. Uses a tissue-to-cell-type graph and a gene-to-cell-type graph, retrieves similar reference cells, then refines with markers.
- **ReCellTy.** arXiv:2505.00017 and bioRxiv 10.1101/2025.04.23.650201 (2025). GraphRAG in Neo4j over a graph built from CellMarker2.0. GPT-4o-mini rose from 0.50 to 0.67 on a human-scored Azimuth set.
- **scHilda.** *PLoS Computational Biology* 22(5):e1014291 (2026), PMC13175463. Knowledge graph from CL, GO, Reactome and MSigDB. It beat GPTCelltype by 0.049 and CellTypeAgent by 0.038 in agreement score.
- **GATHER.** SIGIR 2026, doi:10.1145/3805712.3809935, arXiv:2605.06403. Walks a 20-database knowledge graph from all input genes at once, weights genes IDF-style, and returns the nodes many genes converge on as evidence.
  - Exact match was 27.45% on Immune (34 fine types) and 59.64% on Tabula Sapiens lung. It used one GPT-4o-mini call per cell and beat the graph-reasoning RAG baselines ToG, ToG-2, RoG and PoG.
- **CellTypeAI.** *Bioinformatics* (2026), PMC13501282. Local models with a hand-built JSON marker file used as RAG.
- **scMarkerAgent.** Research Square rs-9224570 (2026). 890k marker statements mined from 294k PMC full texts, each with a supporting sentence.
- **MAT-Cell.** arXiv:2604.06269 (2026). Adaptive RAG plus verifier agents that debate over reasoning trees.
- **RAGCell.** arXiv:2609.14147, posted 12 Sep 2026. Uses retrieved LLM text embeddings as a training signal for a small encoder.

**Foundation models and embedding methods**
- **scGPT and Geneformer.** Kedzierska et al., *Genome Biology* (2025), doi:10.1186/s13059-025-03574-x. Used zero-shot, their embeddings often separated cell types worse than plain HVGs, scVI or Harmony.
- **GenePT.** Chen and Zou, *Nature Biomedical Engineering* (2025), PubMed 39643729. Embeds NCBI gene summaries with GPT-3.5 and weights them by expression. Often matches or beats Geneformer.
- **C2S-Scale.** bioRxiv 10.1101/2025.04.14.648850 (2025). Models from 410M to 27B parameters, with open weights on Hugging Face (vandijklab/C2S-Scale-Gemma-2-27B). It can caption clusters.
- **Cell-o1.** *Bioinformatics* 42(5):btag208 (2026), arXiv:2506.02911. Qwen2.5-7B trained with SFT and then RL. On its batch-level puzzle task, o1 reached only 19%.
- **scAgent.** arXiv:2504.04698 (2025). DeepSeek-R1 plans, scGPT embeds, and a vector database retrieves similar cells. Macro-F1 89.31%.

**Benchmarks**
- **SOAR.** *Briefings in Bioinformatics* 26(6):bbaf622 (2025), arXiv:2412.02915.
- **CellVerse.** NeurIPS 2025 Datasets and Benchmarks, arXiv:2505.07865. 14 LLMs from 160M to 671B parameters. Overall weak results, with drug response near random.
- **DeepCellSeek.** Covered above.

## Q1. Does any method use regulatory network or driver/TF activity as the query or evidence?

**No method I found does this.** Every annotation method above sends marker genes or DEGs to the LLM or the retriever. A few add other signals:
- **Pathway enrichment:** GPTBioInsightor and CyteType.
- **Full DE statistics:** CASSIA Annotation Boost.
- **Filtered marker lists:** LLM-scCurator (bioRxiv 10.64898/2025.12.28.696778) and CELLIA (*BioData Mining* 19:76, 2026, doi:10.1186/s13040-026-00593-7).

The nearest cases:
- **scPilot** gives the LLM pySCENIC TF-gene pairs, but only for a separate task: deciding whether a TF regulates a gene. Its annotation task still uses the top 10 markers.
- **CellMaster's** rationales sometimes cite regulators (for example, IRF8 for pDC). That comes from the LLM's own knowledge, not from a computed activity score.
- **GeneAgent** includes Transfac among its verification databases, but it explains gene sets, not clusters.
- **EpiAgent** (*Nature Methods* 22:2316, 2025, doi:10.1038/s41592-025-02822-z) is a scATAC foundation model, not an LLM agent.

I found nothing that feeds SCENIC AUCell, VIPER/metaVIPER, DoRothEA/decoupleR, or NetBID/scMINER driver activity into an LLM or retriever as evidence. That includes "hidden drivers" with little differential expression. This is exactly the output scMINER produces.

## Q2. Does any method feed the LLM verdict back into the clustering?

Yes, a few do. This is the most crowded part of your planned contribution.
- **GPTAnno.** Song, ..., Hou, bioRxiv 10.1101/2025.11.27.690951 (2025), PMC12687786. It is from the GPTCelltype author's lab.
  - Runs clustering at resolutions 0.1 to 0.7 in steps of 0.2 and queries GPT-5 repeatedly for each cluster.
  - Maps the answers to CL and scores each resolution on consistency, robustness and dispersion on the ontology tree. It picks the resolution with the best combined score.
  - It favours coarse resolutions, then subclusters hierarchically. No open-weight models were tested.
- **HypoGeneAgent.** Genentech, arXiv:2509.09740 (2025).
  - Leiden at resolutions 0.1 to 1.0 on K562 CRISPRi Perturb-seq data.
  - GPT-o3 retrieves from GO, KEGG and PubMed and proposes up to 5 GO hypotheses for each cluster.
  - The resolution score combines agreement among a cluster's hypotheses (ICS, weight 1/3) with how different the top hypothesis is from those of other clusters (1 − ICD, weight 2/3), using text embeddings.
  - It calls itself "the first framework that closes the loop." It targets GO programs, not cell types, and the paper calls itself a preliminary test.
- **CellMaster.** arXiv:2602.13346 (2026).
  - Each loop outputs a label, a confidence score and a recommended next step (split/merge, finer markers, zoom-in). Zoom-in reclusters the chosen clusters at resolution 0.1 to 2.0.
  - A human accepts the suggestions. Average score was 0.602 across 9 datasets, +0.071 over the best baseline. Human feedback added +0.115.
- **Tools that only flag clusters:** mLLMCelltype, LICT, CyteType (heterogeneity) and CASSIA (its subclustering agent annotates subclusters the user supplies).
- **CLUEY.** bioRxiv 10.1101/2024.11.14.623697 (2024). It picks the number of clusters recursively by correlating clusters with atlas cell-type profiles. It uses no LLM, but it is a natural baseline.

I found no work that:
- merges or splits pairs of clusters based on retrieval-grounded label equivalence,
- works with mutual-information clustering, or
- scores the resulting partitions with ARI/NMI against over-clustered baselines. I could not confirm whether GPTAnno reports ARI.

## Q3. Does any method use retrieval to verify or cite labels, and how is faithfulness measured?

- **GeneAgent** has the most rigorous check. Of 15,903 claims, 84% were supported and 8% refuted. Experts reviewed 132 claims from 10 gene sets and found 92% of the verification decisions correct. The database a gene set came from is masked during verification to avoid leakage.
- **CellTypeAgent** grounds labels in CELLxGENE expression, but only final accuracy is reported.
- **CASSIA** quality scores go with correctness (Kruskal-Wallis p<0.01). Agreement across repeated runs also flags errors.
- **LICT** checks that proposed markers are actually expressed in the data.
- **CyteType** reports that confidence goes with similarity (F=23.88). UNVERIFIED detail.
- **scHilda** shows traceable evidence chains and admits explanations sometimes do not match answers. It gives no metric for this.
- **ReCellTy** shows the retrieved entities and measures no faithfulness.
- **GATHER** counts evidence items and lists expert validation of retrieved nodes as future work.
- **CellTypeAI** reports that RAG added 10 to 15% accuracy, with no citations.
- **scMarkerAgent** attaches a supporting sentence to each database entry, not to each annotation.

**Gap:** no cell-annotation paper measures citation precision/recall, NLI-based support of each sentence in its rationale, or ALCE-style attribution against retrieved PubMed passages.

## Q4. Datasets, ontologies, metrics and typical accuracy

**Datasets**
- The GPTCelltype collection (Azimuth, GTEx, HCL, MCA, Tabula Sapiens, BCL, colon and lung cancer) is reused by CASSIA, GPTAnno, scHilda and GPTBioInsightor.
- Other common choices: PBMC3k, Tabula Sapiens v2, CELLxGENE, HypoMap.

**Metrics**
- The dominant one is CL-based scoring of full 1, partial 0.5, mismatch 0 (GPTCelltype, CASSIA, mLLMCelltype, CellTypeAgent, CellMaster). GPTAnno also counts a child term as 1.
- Others:
  - an LLM judge (AnnDictionary)
  - embedding similarity (CyteType, ReCellTy, HypoGeneAgent)
  - BLEU and F1 (SOAR)
  - exact match plus CL ancestor match (GATHER)
  - a custom 3-point manual rubric (CellTypeAI)

**Typical accuracy**

| Setting | Result |
|---|---|
| Broad types, top API models | 80 to 90% |
| mLLMCelltype vs GPTCelltype | 77.3% vs 61.3% |
| scPilot with o1 | 0.79 (PBMC3k), 0.52 (liver) |
| CellMaster | 0.60 average |
| Fine immune subtypes, cell level (GATHER) | 27% exact |

Accuracy drops for disease states, low-heterogeneity tissue (LICT: 43.8 to 48.5%) and subtypes.

## Q5. Open-weight models, and whether a local 30B to 70B model can compete

**Evidence by paper**
- **SOAR (2025, older models, no RAG).** Clear gap. Zero-shot chain-of-thought F1: GPT-4o 57.36, Qwen2-72B 46.38. Llama-3-70B was the weakest general model. Mixtral-8x22B had the best open BLEU (28.19).
- **CASSIA preprint (100 cell types).** Claude 3.5 Sonnet 0.92, GPT-4o 0.88, Llama-3.2-90B 0.82.
- **AnnDictionary.** Llama 3.1 405B was near the top, around 82%. The 70B and 8B models fell in lower tiers (exact values UNVERIFIED).
- **scHilda.** DeepSeek-V3.2 came close to o3.
- **CyteType.** DeepSeek R1, Qwen3 and Kimi K2 were modestly but significantly below closed models.
- **CellTypeAgent.** DeepSeek-R1 gained 5.1% with database verification.
- **CellTypeAI (with RAG).**
  - PBMC3k: Qwen3-32B 87.5% and Qwen3-235B 94.4%, against Gemini 2.5 Pro at 86.6%.
  - Stroke PBMC: Qwen3-32B 86.1% against Gemini at 84.0%.
  - These are small, custom-rubric evaluations.
- **scPilot.** Gemma-3-27B underperformed and was dropped.
- **Not open-weight:** Qwen2.5-Max in mLLMCelltype is API-only.

**Hardware estimate (my own calculation, not from a paper)**
- Qwen3-32B in BF16 needs about 65 GB and fits across your 2x48 GB with tensor parallelism.
- A 70B model needs FP8 (about 70 GB, tight on KV cache) or 4-bit AWQ (about 40 GB).
- gpt-oss-120b in MXFP4 (about 61 GB) should fit. I found no annotation benchmark for it.
- Qwen3-235B does not fit in VRAM. CPU offload with your 754 GB RAM would work but be slow.

**Verdict:** a local 32B to 70B model with RAG or agent scaffolding is plausibly within about 5 to 10 points of GPT-4-class models on broad types. It is weaker on fine subtypes. Include a GPT-class API baseline on a subset.

## Gap analysis

1. **Regulator activity as evidence.** Nobody uses computed TF, driver or regulon activity (SJARACNe/NetBID2 as in scMINER, SCENIC, VIPER, DoRothEA) as the query for LLM or RAG annotation. Hidden drivers without DE are invisible to every marker-based method.
2. **Retrieval-verified feedback into clustering.** GPTAnno and HypoGeneAgent select a resolution, and CellMaster suggests split/merge to a human. None of them:
   - makes automatic merge/split decisions for pairs of clusters from retrieval-grounded identity evidence,
   - works on MI-based clustering, or
   - reports partition quality (ARI/NMI) against plain modularity or silhouette selection.
3. **Faithfulness metrics.** No claim-level attribution or citation precision/recall exists for cell annotation. The only measured verifier is GeneAgent's, and it covers gene sets. The IR evaluation methods you already know (attribution, nugget-based scoring) have not been applied here.
4. **Open-weight baselines.** No systematic comparison of local 30B to 70B models with RAG across the standard GPTCelltype collection. The evidence is scattered and small.
5. **Retrieval quality.** Only GATHER, an IR venue paper, treats the many-gene query as a retrieval problem. Dense or ColBERT-style retrieval over PubMed and scMarkerAgent evidence sentences for cluster queries is untested.

**Nearest prior art you must cite:** GPTAnno, HypoGeneAgent, CellMaster, CASSIA, GeneAgent and GATHER.

**Novelty risk is high.** RAGCell was posted on 12 September 2026. A defensible paper would combine items 1, 2 and 3 on scMINER, using local models.

**Queries tried (all returned nothing for Q1 or partial hits for Q2)**
- "LLM cell type annotation transcription factor activity regulon SCENIC input"
- "LLM agent GRN interpretation regulon SCENIC GPT-4 annotation 2025"
- "VIPER OR DoRothEA TF activity LLM cell identity annotation"
- "LLM hidden driver NetBID scMINER"
- "decoupler TF activity LLM prompt annotation"
- "master regulator LLM VIPER cluster annotation"
- "scATAC LLM chromVAR motif annotation"
- "SCENIC ChatGPT regulon cell type"
- "LLM-guided clustering resolution merge split over-clustering"
- "clustering resolution LLM agent subclustering"
- "benchmark hallucination grounding LLM cell type annotation"
- "hallucinated marker genes fabricated citations single-cell"

**UNVERIFIED**
- LICT's PBMC improvement from 21.5% to 90.3% (seen only in a search snippet).
- CyteType's exact numbers and databases.
- CELLIA's Leiden-resolution sensitivity details.
- scLLM-DSC (arXiv:2606.13007): uses LLM semantic anchors inside deep clustering. Details not read.
- Whether scExtract takes the cluster count from the paper.
- The CASSIA README's current model names, and the sc-llm-benchmark.com leaderboard (the domain did not resolve).
- DeepSeq's models (ICML 2025 FM4LS workshop, arXiv:2506.13817).
- Biomni's annotation accuracy.
- Whether CellTypeAgent has since appeared in a peer-reviewed venue.

Local PDFs I read (auto-saved by the fetch tool, not written by me):
- /home/ugdf8/.claude/projects/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/tool-results/webfetch-1789576686097-0e1uet.pdf (CellMaster)
- /home/ugdf8/.claude/projects/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/tool-results/webfetch-1789576683432-i1erbg.pdf (GATHER)
- /home/ugdf8/.claude/projects/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/tool-results/webfetch-1789576774446-dl7q83.pdf (HypoGeneAgent)

Three claude.ai connectors (Gmail, Google Calendar, Google Drive) need authorization in the user's claude.ai connector settings. They were not used here and stay unavailable until then.