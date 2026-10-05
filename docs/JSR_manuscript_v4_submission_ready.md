# Sleep deprivation is associated with cell-type-heterogeneous suppression of intestinal oxidative-phosphorylation transcriptional programs

**Running title:** Sleep loss and intestinal respiration

**Authors:** Qiyun Hu¹, Ziqi Wang¹, Weina Wang¹, Hu Zhang¹*

**Affiliation:** ¹The Second Affiliated Hospital of Zhejiang Chinese Medical University (Xinhua Hospital of Zhejiang Province), Hangzhou, Zhejiang, China

**Corresponding author:** Hu Zhang, The Second Affiliated Hospital of Zhejiang Chinese Medical University (Xinhua Hospital of Zhejiang Province), 318 Chaowang Road, Gongshu District, Hangzhou 310005, Zhejiang Province, China. Email: 20234028@zcmu.edu.cn. ORCID: None / not available.

## Abstract

Sleep loss has been linked to intestinal injury and impaired epithelial repair, but the cellular distribution and cross-study reproducibility of the associated transcriptional programs remain uncertain. To reduce dataset-selection bias, we systematically searched public transcriptomic resources and identified three eligible intestinal sleep-deprivation studies. GSE289089 served as a single-cell discovery dataset, with sample-level pseudobulk analysis used to resolve epithelial heterogeneity. Negative oxidative-phosphorylation enrichment was strongest in intestinal stem cells, transit-amplifying progenitors, tuft cells, goblet cells, and enterocytes, whereas several other epithelial populations showed weaker or nonsignificant results. PRJNA1167170, an independent proximal-colon bulk RNA-sequencing dataset, showed strong negative enrichment of the same program. A second independent colon dataset, PRJNA862187, had a negative oxidative-phosphorylation estimate but did not show statistically significant enrichment (normalized enrichment score −0.602; nominal P = 1.000; false-discovery rate = 1.000). Among five prespecified secondary mitochondrial pathways, mitochondrial translation showed statistically supported negative enrichment in GSE289089 intestinal-stem-cell pseudobulk and in both bulk datasets. Gene-level concordance was strongest between GSE289089 and PRJNA1167170 and was weak for comparisons involving PRJNA862187. Thus, a shared negative direction did not constitute uniform statistical or gene-level replication. These secondary transcriptomic analyses provide no direct mitochondrial functional measurement and do not establish causality.

## Keywords

crypt homeostasis; cellular bioenergetics; biological-sample inference; ranked-list analysis; murine models; public-data reuse

## Introduction

Sleep disruption affects physiology beyond the nervous system, and the intestine has emerged as a potentially important peripheral target. In experimental animals, severe sleep loss can produce gut oxidative stress; in flies, preventing intestinal reactive-oxygen-species accumulation can rescue lethality associated with profound sleep deprivation, while mice also accumulate gut reactive oxygen species after sleep loss (Vaccaro et al., 2020). Mouse studies have additionally associated sleep deprivation with colonic mucosal injury, inflammatory changes, and microbiota alterations (Gao et al., 2019; Park et al., 2020). These findings are paradigm- and species-dependent. Five hours of acute sleep deprivation caused only subtle microbiota changes in mice (El Aidy et al., 2020), and a randomized crossover study in healthy young men found reduced microbial richness but no change in measured intestinal permeability after three nights of severe sleep restriction (Karl et al., 2023). Intestinal responses to sleep loss should therefore not be treated as uniform across exposure models or extrapolated directly from rodents to humans.

Recent animal studies have begun to define specific epithelial mechanisms. Sleep deprivation reduced goblet-cell abundance and mucin expression and activated an endoplasmic-reticulum-stress response in mice (Li et al., 2024). A separate study linked chronic sleep deprivation to colonic metabolic disruption, alpha-ketoglutarate accumulation, prolyl hydroxylase 2-dependent degradation of hypoxia-inducible factor 1 alpha, and impaired repair (Zhang et al., 2024). More recently, acute sleep deprivation was reported to activate a vagal acetylcholine–enterochromaffin serotonin circuit, with serotonin receptor 4-dependent oxidative stress in intestinal stem cells (Zhang et al., 2026). These studies establish biologically specific sleep–gut pathways, but their mechanisms differ and do not show that one epithelial population or mechanism accounts for all intestinal responses.

The intestinal epithelium is well suited to resolving heterogeneity because absorptive, secretory, chemosensory, and continuously renewing compartments coexist within the same tissue. Whole-tissue profiles can identify a pathway direction but cannot determine whether it originates broadly or from a restricted population. Conversely, single-cell data provide cellular resolution but yield misleading precision when cells rather than animals are treated as replicates. Combining cell-resolved discovery, sample-level inference, and independent tissue-level validation can test reproducibility without implying that whole tissue reproduces the cellular origin of a signal.

Public transcriptomic studies can provide such complementary resolution, but selecting datasets after inspecting their results can exaggerate apparent agreement. Repository records also commonly duplicate one experiment across GEO, SRA, BioProject, ENA, publications, and preprints, while incomplete accessions may appear usable from a title alone. A study-level search with explicit eligibility and exclusion rules is therefore important for distinguishing reproducibility from selective dataset inclusion. It also creates a record of eligible studies whose designs differ from the discovery experiment and may expose genuine context dependence rather than technical failure.

Mitochondrial oxidative phosphorylation provides a plausible link between epithelial energy metabolism and renewal. In the small-intestinal crypt, Lgr5-positive crypt-base-columnar cells show high mitochondrial activity, and Paneth-cell-derived lactate sustains their enhanced oxidative phosphorylation and stem-cell function (Rodríguez-Colman et al., 2017). Nevertheless, transcriptional enrichment cannot establish organelle dysfunction, and study selection can bias apparent reproducibility. We therefore sought to (1) characterize epithelial cell-type heterogeneity, (2) test independent transcriptomic reproducibility, (3) reduce dataset-selection bias through systematic identification of eligible public studies, and (4) determine whether mitochondrial transcriptional direction persisted across heterogeneous experimental contexts. The study adds systematic public-dataset identification, biological-sample-level pseudobulk resolution, cross-study pathway validation, gene-level concordance where supported, and explicit reporting of context dependence.

## Methods

### Study design and analysis sequence

This study was a secondary reanalysis of publicly available transcriptomic data from experimental mouse studies. GSE289089/PRJNA1221288 was the single-cell discovery dataset and contained two control and two 48-hour sleep-deprivation samples. PRJNA1167170 was the first independent validation dataset and comprised four control and four sleep-deprivation proximal-colon bulk RNA-sequencing samples. PRJNA862187 was identified as an additional eligible study during a subsequent systematic public-data search and was analysed as a second independent colon validation dataset. The datasets were modelled separately and were not pooled, integrated, or treated as removable technical batches. Biological samples, not individual cells, were the units of replication.

The GSE289089 discovery analysis and PRJNA1167170 validation were completed before the systematic search identified PRJNA862187. Before downloading or examining PRJNA862187 expression results, its analysis plan, factorial contrasts, primary pathway, expected direction, secondary pathways, and interpretation rules were frozen in `PRJNA862187_analysis_plan_frozen.md`. The three datasets were therefore not simultaneously preregistered. Following its identification in the discovery dataset, `HALLMARK_OXIDATIVE_PHOSPHORYLATION` was designated, before analysis of each independent validation dataset, as the primary validation pathway.

### Public dataset identification and selection

To reduce dataset-selection bias, we systematically searched NCBI GEO, SRA, BioProject, PubMed, and Europe PMC through 4 October 2026. ENA and BioStudies/ArrayExpress were used to cross-check accession status and public availability. Searches combined sleep-exposure concepts (sleep deprivation, restriction, loss, disruption, and disturbance), intestinal concepts (intestine, gut, colon, small intestine, epithelium, crypt, intestinal stem cell, goblet, and enterochromaffin), and transcriptomic concepts (RNA-seq, transcriptome, scRNA-seq, single-nucleus RNA sequencing, and microarray). Repository accessions, publications, and preprints describing the same experiment were merged into one study-level record.

Predefined eligibility criteria required an explicit sleep-deprivation or sleep-restriction exposure, a contemporaneous control group, intestinal or colonic tissue or an intestinal epithelial population, transcriptomic profiling, at least two independent biological samples per condition, interpretable group labels, accessible public data at the cutoff, and sufficient information to construct a sleep-deprivation-versus-control comparison. Standardized exclusions covered non-intestinal tissue, absence of a sleep contrast or control, non-transcriptomic assays, insufficient replication, inaccessible data, unresolved group labels, internally inconsistent library metadata, duplicate studies, and absence of reusable expression data. Unknown fields were not inferred. This was a targeted public-dataset identification and selection audit, not a PRISMA systematic review.

The complete search log retained exact query wording, overlapping database hit counts, accession checks, and decision notes. Nineteen exclusions comprised seven non-intestinal studies, two without an eligible sleep-deprivation or sleep-restriction contrast, three with non-transcriptomic assays only, five without accessible public data, one with unresolved group labels, and one without reusable expression data. Duplicate repository records were merged rather than counted as exclusions. `POTENTIALLY_ELIGIBLE` was reserved for a biologically relevant design whose repository metadata did not permit confident expression analysis.

### Single-cell preprocessing, annotation, and pseudobulk modelling

Four GSE289089 count matrices were used to construct Seurat objects (Stuart et al., 2019). The contemporaneous record specified `CreateSeuratObject(min.cells = 3, min.features = 200)`. All samples used the same quality-control thresholds: `nFeature_RNA` 500–6,500 and `percent.mt < 25`, with no hard `nCount_RNA` cutoff. Quality control reduced 36,043 cells to 30,580. Doublets were identified using the recovered call `scDblFinder(sce, samples = "sample")` (Germain et al., 2022). No additional arguments were present in the recovered call. Because the original package version was unavailable, the exact effective default values could not be retrospectively reconstructed. Removal of 2,909 classified doublets left 27,671 singlets.

The RNA assay was normalized with `LogNormalize` and a scale factor of 10,000; 2,000 variable features were selected using the variance-stabilizing-transformation method. Data were centered and scaled without a recorded regression variable. Principal-component analysis used 50 components and seed 42. The neighbor graph used components 1–20, `k.param = 20`, Annoy Euclidean search, and 50 trees. Louvain clustering used algorithm 1, resolution 0.5, 10 starts, 10 iterations, and seed 123. UMAP used components 1–20, the `uwot` implementation, cosine distance, 30 neighbors, minimum distance 0.3, spread 1, and seed 123.

Frozen identities for clusters 0–17 were retained. The original marker call was `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)`. Because its saved output was unavailable, the same call and canonical panel were applied once to the frozen object solely to reconstruct annotation evidence. This did not alter normalization, clustering, identities, or downstream results. The original cell-cycle record used Seurat `cc.genes.updated.2019`; 43 S-phase and 54 G2M genes were matched after conversion to mouse symbols. `CellCycleScoring` used `set.ident = FALSE`, and cell-cycle scores were not included as covariates or inferential outcomes.

Clusters were mapped to intestinal stem cells; transit-amplifying/cycling progenitors; goblet/secretory cells; Paneth lineage; enterocytes; enteroendocrine cells; Tph1-positive serotonergic enterochromaffin cells; and tuft cells. Stress-responsive, interferon-responsive/regenerative, and enteroendocrine-progenitor states remained exploratory. Cluster 17 contained 34 immune cells and was excluded from epithelial comparisons. Raw counts were aggregated by biological sample within cell group using `AggregateExpression`; formal inference required at least 20 cells in every sample.

Differential expression used edgeR (Robinson et al., 2010). Within each eligible group, `DGEList`, `filterByExpr`, `calcNormFactors`, `estimateDisp`, `glmQLFit`, and `glmQLFTest` were applied with model `~ condition`. Source symbols were harmonized to GENCODE Mouse M39 using stable Ensembl identifiers. Duplicate source feature symbols in the 10x input were deterministically disambiguated with `make.unique` so that each Seurat row remained mapped to one Ensembl identifier; counts and edgeR models remained indexed by Ensembl identifier. For symbol-indexed enrichment and concordance, multiple Ensembl identifiers mapping to the same current symbol were reduced by retaining the row with the highest edgeR log counts per million; values were not summed or averaged. Benjamini–Hochberg FDR values were reported. This approach preserved animal-level replication and did not estimate differential cell abundance.

Aggregation and modelling were repeated independently for each eligible cell group, allowing contributing cell numbers and library sizes to differ among groups and samples while retaining the animal as the replicate. Sample-level cell counts were used only for eligibility checks, descriptive displays, and the frozen downsampling analysis. The stress-responsive, interferon-responsive/regenerative, and enteroendocrine-progenitor states were not silently merged into lineage results and were reported only as exploratory populations.

### Gene-set enrichment and sensitivity analyses

Native-mouse MSigDB 2026.1.Mm sets were obtained with msigdbr 26.1.1 (Liberzon et al., 2015). `HALLMARK_OXIDATIVE_PHOSPHORYLATION` was the primary validation pathway after discovery. Respiratory electron transport, mitochondrial translation, ATP formation by chemiosmotic coupling, the tricarboxylic-acid cycle, and mitochondrial fatty-acid beta oxidation were prespecified secondary validation pathways. Gene-set enrichment used `fgseaMultilevel` (Korotkevich et al., 2021) with `eps = 0`, minimum size 10, maximum size 500, and edgeR log fold change as the primary ranking statistic (Subramanian et al., 2005). Because fgsea evaluates enrichment conditional on the fitted ranked gene list, pathway-level FDR values were not interpreted as a substitute for biological replication.

A frozen sensitivity ranking used `sign(logFC) × sqrt(F)`. Six cell groups were also downsampled within each sample to the smallest sample-specific cell count for 20 iterations using master seed 20261003. Regeneration-associated pathways were exploratory and retained with their mixed directions. These analyses tested ranking and cell-number sensitivity but did not add biological replicates.

### First independent proximal-colon validation dataset

PRJNA1167170 comprised paired-end reads from four control and four sleep-deprivation proximal-colon tissue samples. FastQC and MultiQC were used for quality control; no sample was excluded. Salmon 2.7.0 (Patro et al., 2017) quantified transcripts against a GENCODE Mouse M39/GRCm39 index with automatic library detection. All libraries were inferred as IU. Transcript estimates were imported with tximport 1.40.0 using the corresponding M39 transcript-to-gene map and `countsFromAbundance = "no"` (Soneson et al., 2015). Gene-level edgeR modelling used `~ condition` and the quasi-likelihood workflow. The primary OXPHOS validation used log-fold-change ranking and `fgseaMultilevel`.

### Second independent validation dataset

PRJNA862187 contained 12 independent colon RNA-sequencing samples in a 2×2 sleep-by-diet design: exercise control with standard chow (EC_SCD), sleep deprivation with standard chow (SD_SCD), exercise control with high-fat diet (EC_HFD), and sleep deprivation with high-fat diet (SD_HFD), with three biological replicates per group (Lee et al., 2023). In the source study, EC denoted exercise control, not an undisturbed cage control: EC mice received the same total walking distance using longer wheel on/off intervals that allowed longer uninterrupted sleep, whereas the sleep-deprivation schedule used repeated brief wheel movements to interrupt sleep. The EC condition was intended to control movement-related nonsleep effects of the forced-wheel paradigm. Paired-end 151-bp reads were quantified with Salmon 2.7.0 against the same M39/GRCm39 reference using automatic library detection. All libraries were inferred as ISR, and mapping rates ranged from 88.34% to 93.05%. No sample was excluded. Transcript estimates were imported through the same framework with `countsFromAbundance = "no"`.

The edgeR design was `~0 + group`. The primary contrast estimated the average sleep-deprivation effect across diets: `0.5×SD_SCD + 0.5×SD_HFD − 0.5×EC_SCD − 0.5×EC_HFD`. Secondary contrasts were `SD_SCD − EC_SCD` and `SD_HFD − EC_HFD`. The exploratory interaction was `(SD_HFD − EC_HFD) − (SD_SCD − EC_SCD)`. The primary validation endpoint remained Hallmark OXPHOS with a prespecified expected direction of NES below zero. The five mitochondrial pathways listed above were secondary endpoints.

### Cross-dataset concordance

For all native-mouse Hallmark OXPHOS genes available in each fitted dataset, sleep-deprivation log fold changes were compared without significance filtering. Pairwise Spearman correlations were calculated between GSE289089 intestinal-stem-cell pseudobulk, PRJNA1167170 proximal-colon bulk, and the PRJNA862187 average sleep-deprivation contrast. When multiple tested Ensembl identifiers shared one current gene symbol, the identifier with the highest log counts per million was retained before symbol-level matching. We counted genes measured and negative in all three datasets. Pathway NES values were summarized separately because tissue-level validation and intestinal-stem-cell pseudobulk differ in biological resolution, and pathway-level direction does not imply cell-type-matched replication.

### Software and reproducibility

Frozen analyses used R 4.6.1, Seurat 5.5.1, SeuratObject 5.4.0, edgeR 4.10.5, fgsea 1.38.0, msigdbr 26.1.1, tximport 1.40.0, data.table 1.18.6.1, ggplot2 4.0.3, and Salmon 2.7.0. Source tables, provenance records, scripts, and figure-source mappings were retained locally. No frozen model or enrichment result was recomputed during manuscript integration.

## Results

### Systematic public-data search identified three eligible transcriptomic studies

After study-level deduplication, 23 unique candidate studies underwent full metadata assessment. Three were eligible: GSE289089/PRJNA1221288, PRJNA1167170, and PRJNA862187. Nineteen were excluded, and PRJNA1395688 (PMID 41094491) remained potentially eligible because public raw reads were linked to internally inconsistent sequencing-library metadata. PRJNA1019685 and an additional 2026 intestinal transcriptomic study were excluded because no reusable public expression data were available at the search cutoff. The complete registry, study information, and exclusion reasons are provided in the Supporting Information.

GSE289089 was used for single-cell discovery, PRJNA1167170 for the first independent proximal-colon validation, and PRJNA862187 for a second independent colon validation under the analysis plan frozen before data download. This sequence preserved the distinction between discovery, validation, and subsequent systematic identification.

### Single-cell analysis showed epithelial heterogeneity

GSE289089 contained two control and two sleep-deprivation biological samples. Quality control and doublet removal retained 27,671 singlets in 18 frozen clusters (Figure 1). Reconstructed marker evidence supported the frozen annotations without relabelling. Sample-level cell-group proportions were descriptive only.

The primary analysis showed significant negative OXPHOS enrichment in intestinal stem cells (NES = −2.032, FDR = 2.14 × 10^−7), transit-amplifying progenitors (NES = −1.784, FDR = 4.33 × 10^−6), tuft cells (NES = −2.271, FDR = 2.89 × 10^−11), goblet cells (NES = −1.507, FDR = 0.00619), and enterocytes (NES = −1.476, FDR = 0.00276) (Figure 2). Paneth-lineage cells had weaker nonsignificant negative enrichment (NES = −1.103, FDR = 0.252). Enteroendocrine cells (NES = −0.497, FDR = 1.000) and serotonergic enterochromaffin cells (NES = −0.618, FDR = 1.000) also lacked significant primary enrichment. The response was therefore heterogeneous and not intestinal-stem-cell-specific.

The discovery dataset contained only two biological samples per condition, and contributing cell counts were imbalanced in some populations. The most marked example was the goblet group, with 77 cells in 2D_1 and 947 in 2D_2. These results are sample-level transcriptional comparisons, not estimates of differential cell abundance.

Reconstructed canonical-marker evidence supported the frozen lineage and state labels without prompting relabelling. Stem-associated, cycling, Paneth, goblet, enterocyte, enteroendocrine, tuft, and immune programs were represented in their corresponding clusters. Absence of an individual marker at the original threshold was not used to replace a working identity when the broader pattern remained supportive. Because the original marker output was unavailable, these displays document annotation provenance rather than a new classification analysis.

### The first independent proximal-colon dataset provided strong validation

PRJNA1167170 contained four control and four sleep-deprivation proximal-colon samples. Mapping rates ranged from 89.902% to 95.196%, and all libraries were inferred as IU. Hallmark OXPHOS showed strong negative enrichment (NES = −2.569, FDR = 6.86 × 10^−15; Figure 3). All five prespecified secondary pathways were also significantly negative: TCA cycle (NES = −1.635, FDR = 0.0446), ATP formation by chemiosmotic coupling (NES = −2.189, FDR = 2.16 × 10^−5), mitochondrial fatty-acid beta oxidation (NES = −1.546, FDR = 0.0478), mitochondrial translation (NES = −2.397, FDR = 3.24 × 10^−8), and respiratory electron transport (NES = −2.698, FDR = 9.31 × 10^−14). This was independent proximal-colon tissue validation, not cell-type-matched replication.

### A second eligible colon dataset showed a negative but statistically unsupported OXPHOS estimate

All 12 PRJNA862187 samples completed quantification, with no technical outlier excluded. For the primary average sleep-deprivation effect across diets, OXPHOS enrichment was negative but not statistically significant (NES = −0.602, P = 1.000, FDR = 1.000; Figure 3). Diet-specific estimates were also negative but nonsignificant under standard chow (NES = −0.401, FDR = 1.000) and high-fat diet (NES = −0.789, FDR = 0.955). PRJNA862187 therefore had the same estimated direction as the earlier datasets but did not provide statistical replication of OXPHOS enrichment.

Among the five prespecified secondary pathways, mitochondrial translation showed the clearest PRJNA862187 signal. It was negatively enriched for the average sleep-deprivation effect (NES = −1.748, FDR = 0.00125) and under high-fat diet (NES = −2.255, FDR = 2.15 × 10^−7), whereas the standard-chow estimate was positive and nonsignificant (NES = 0.887, FDR = 0.987). Mitochondrial translation remained a secondary endpoint. Respiratory electron transport was negative under high-fat diet (NES = −1.407, FDR = 0.0491), while the other average-effect secondary results were not statistically clear.

No gene-level sleep-by-diet interaction survived FDR correction among 17,900 tested genes. The OXPHOS interaction was nonsignificant (NES = 0.382, FDR = 1.000), providing no clear evidence of OXPHOS interaction. An exploratory pathway-level interaction was observed for mitochondrial translation (NES = −1.765, FDR = 0.000733), but no gene-level interaction survived FDR correction.

### Pathway direction was more consistent than gene-level concordance across three datasets

Among Hallmark OXPHOS genes, GSE289089 and PRJNA1167170 showed positive concordance (Spearman ρ = 0.337, P = 1.86 × 10^−6). Concordance was not clear between GSE289089 and PRJNA862187 (ρ = 0.079, P = 0.276) and was weakly positive between PRJNA1167170 and PRJNA862187 (ρ = 0.154, P = 0.0321) (Figure 4). Of 195 Hallmark genes, 191 were measurable in all three datasets and 82/191 (42.93%) were negative in all three. Thus, the negative pathway direction was more consistent than the magnitude or ranking of individual-gene effects.

For the originally completed two-dataset comparison, 138/191 genes (72.25%) were negative in both GSE289089 intestinal-stem-cell and PRJNA1167170 proximal-colon results. Their formal leading edges contained 118 and 133 genes, respectively; 90 overlapped within a 161-gene union (Jaccard index 0.559). These results support concordance between discovery and the first validation dataset but do not extend strong gene-level replication to PRJNA862187.

At pathway level, OXPHOS, respiratory electron transport, mitochondrial translation, ATP formation by chemiosmosis, and the TCA cycle were negative in GSE289089 intestinal-stem-cell pseudobulk and in both bulk datasets, although most PRJNA862187 estimates were weak. Mitochondrial translation was the only prespecified secondary pathway with statistically supported negative enrichment in the GSE289089 intestinal-stem-cell analysis, PRJNA1167170 proximal-colon tissue, and PRJNA862187 colon tissue. Mitochondrial fatty-acid beta oxidation was negative in the first two datasets but positive and nonsignificant in PRJNA862187.

### Sensitivity analyses retained lineage-specific uncertainty

The alternate signed-statistic ranking reproduced significant negative OXPHOS enrichment in intestinal stem cells, transit-amplifying progenitors, tuft cells, goblet cells, and enterocytes and additionally yielded negative enrichment in Paneth-lineage cells. Enteroendocrine cells remained nonsignificant and changed direction; serotonergic enterochromaffin cells remained negative but nonsignificant. Across 20 equal-cell downsampling iterations, the five strongest groups remained negative with FDR below 0.05 in every iteration under both rankings. Regeneration-associated programs remained mixed rather than showing coherent suppression. These sensitivity analyses did not add animals or constitute independent replication.

Paneth-lineage enrichment under the alternate ranking was negative (NES = −1.453, FDR = 0.00712), whereas the enteroendocrine result was positive and nonsignificant (NES = 0.977, FDR = 0.523). Serotonergic enterochromaffin enrichment remained negative but nonsignificant (NES = −1.037, FDR = 0.362). During downsampling, the serotonergic enterochromaffin result was negative in 90% of log-fold-change-ranked and 95% of signed-statistic-ranked iterations but reached FDR below 0.05 in only 0% and 10%, respectively. This distinction prevents frequent negative direction from being presented as stable statistical support.

## Discussion

Systematic searching identified three eligible public intestinal sleep-deprivation transcriptomic studies. The strongest evidence for a negative OXPHOS transcriptional response came from GSE289089 and PRJNA1167170. PRJNA862187 had the same negative primary estimate but no statistical support for enrichment. The overall finding is therefore a shared direction with heterogeneous statistical evidence, not uniformly strong replication across contexts. Cell-resolved discovery further showed that the response was prominent in intestinal stem cells, transit-amplifying progenitors, tuft cells, goblet cells, and enterocytes but weaker or nonsignificant in several other populations. The signal is neither uniform across epithelial lineages nor specific to intestinal stem cells.

This pattern complements, but does not reproduce mechanistically, prior animal experiments. Goblet-cell endoplasmic-reticulum stress (Li et al., 2024), alpha-ketoglutarate/PHD2/HIF1alpha-associated repair impairment (Zhang et al., 2024), and a vagal acetylcholine–enterochromaffin serotonin circuit affecting intestinal stem cells (Zhang et al., 2026) provide relevant biological frameworks. The present secondary analysis did not measure these metabolites, proteins, neural signals, or causal pathways. Its contribution is to define the distribution and cross-study direction of transcriptional programs, not to adjudicate among those mechanisms.

PRJNA862187 underscores context dependence. Its factorial diet design, sleep paradigm, colon preparation, and greater within-group expression heterogeneity differed from the discovery and first validation studies. These features may contribute to its weaker OXPHOS signal and weak gene-level concordance, but the present data cannot assign causality to any one difference. The negative OXPHOS estimates under both diets and nonsignificant OXPHOS interaction do not provide clear evidence that diet modified the primary pathway response. Conversely, absence of a significant interaction is not evidence that diet has no influence.

Mitochondrial translation was the most consistent prespecified secondary mitochondrial pathway across the GSE289089 intestinal-stem-cell analysis and the two bulk datasets. This pattern may indicate that regulation of mitochondrial protein synthesis is a cross-study component of the transcriptional response. However, it was a secondary endpoint, and transcriptomic enrichment does not establish impaired mitochondrial translation or protein synthesis. Its exploratory interaction in PRJNA862187 suggests possible dietary modification but requires independent confirmation, particularly because no gene-level interaction survived FDR correction.

The three-dataset comparison also separates pathway replication from gene-level replication. GSE289089 and PRJNA1167170 showed moderate positive gene-level concordance, whereas the comparison involving PRJNA862187 was weak or unclear. Nonetheless, OXPHOS remained negative in all three datasets and 82 of 191 measurable genes were negative throughout. Pathway-level statistics can remain directionally aligned even when individual-gene rankings show weak concordance. Whole-tissue datasets also mix cell types and cannot be considered intestinal-stem-cell-specific replication.

Pathway FDR values must be interpreted within each fitted ranked list. GSE289089 had only two biological samples per condition, and PRJNA862187 produced a negative OXPHOS estimate with nominal P = 1.000 and FDR = 1.000. Direction and statistical support were therefore reported separately; the PRJNA862187 result was not treated as statistical replication.

Several findings prevent an overly simple interpretation. Paneth and enteroendocrine-lineage primary results were nonsignificant, enteroendocrine enrichment changed direction under the alternate ranking, and regeneration-associated pathways were mixed. Mitochondrial fatty-acid beta oxidation was not directionally consistent in the third dataset. These negative and mixed findings argue against relabelling the phenotype as general mitochondrial dysfunction, which would require direct evidence from respiration, ATP production, membrane potential, mitochondrial structure, protein synthesis, or related functional measurements.

The broader sleep–gut literature also cautions against uniform interpretation. Some mouse paradigms link sleep deprivation to barrier injury, inflammation, and microbiota dysbiosis (Gao et al., 2019; Park et al., 2020), whereas brief acute deprivation produced only subtle microbiota changes (El Aidy et al., 2020). In healthy young men, severe short-term restriction reduced microbial richness without altering measured intestinal permeability (Karl et al., 2023). Differences in duration, severity, host, tissue, and outcome do not invalidate the present mouse transcriptional findings, but they limit translation to other sleep exposures and to humans.

The study has important limitations. GSE289089 contained only two biological samples per condition; pseudobulk analysis prevents cell-level pseudoreplication but cannot create biological replication. Cell numbers were uneven across samples, although equal-cell downsampling supported the strongest signals. PRJNA1167170 and PRJNA862187 were bulk colon datasets rather than cell-type-matched validations, and differences in intestinal segment, preparation, composition, sleep exposure, and diet can affect concordance. The public-dataset search was a targeted identification audit rather than a formal systematic review. PRJNA1395688 remained unresolved because of conflicting library metadata, while reusable transcriptomic data for PRJNA1019685 and an additional 2026 intestinal transcriptomic study were unavailable at the cutoff. All analyses were secondary reanalyses of experimental mouse studies, so association cannot establish causation.

Provenance was also incomplete for parts of the original single-cell workflow. The QC thresholds and scDblFinder call were recovered, but the exact mitochondrial-feature pattern, per-sample execution structure, scDblFinder package version, and internally resolved defaults were unavailable. Marker evidence was reconstructed with the original parameters because the saved marker output was absent. These gaps do not change the frozen cell set but should remain transparent.

The independent datasets have complementary limitations. PRJNA1167170 provided strong pathway evidence but only eight proximal-colon samples, and PRJNA862187 contained three colon samples per factorial cell with visually substantial within-group dispersion. No sample was removed because of its PCA position. Differences between intestinal-stem-cell pseudobulk, proximal-colon expression, and colon expression may reflect cell composition as well as cell-intrinsic responses. Consequently, the analysis can support cross-study direction but cannot localize the bulk signals to a particular epithelial population.

Future studies should combine controlled sleep manipulations with matched intestinal segments, larger biological replication, and cell-type-resolved measurements of respiration and mitochondrial protein synthesis. Factorial designs should be independently replicated before assigning diet-dependent effects. Such experiments should retain the epithelial heterogeneity observed here rather than assuming a single cellular target.

## Conclusion

Across systematically identified public datasets, sleep deprivation was associated with negative estimates of intestinal mitochondrial respiratory transcriptional programs, with substantial variation in effect magnitude and statistical support. Single-cell analysis showed marked epithelial heterogeneity rather than an intestinal-stem-cell-specific response. PRJNA862187 had a negative OXPHOS estimate but no statistically significant enrichment, while mitochondrial translation emerged as the most consistent secondary pathway across the GSE289089 intestinal-stem-cell and two bulk analyses. These transcriptomic associations do not demonstrate mitochondrial dysfunction or causality and require functional validation in experimentally matched tissues and cell populations.

## Data Availability Statement

The datasets analysed are publicly available under GSE289089/PRJNA1221288, PRJNA1167170, and PRJNA862187. PRJNA1395688 was identified as potentially eligible but was not analysed because sequencing-library metadata remained internally inconsistent. Analysis scripts and processed outputs are currently retained in the authors' local project archive and have not been deposited in a public repository. Repository information or accurate access conditions will be finalized before submission: [REPOSITORY/DOI OR ACCESS CONDITIONS TO BE PROVIDED].

## Ethics Statement

Ethical approval was not required for this secondary analysis of publicly available archived datasets. No human participants were recruited and no new animal experiments were performed by the authors of this analysis. Experiment-specific ethical approvals remain the responsibility of the original studies and should be verified from those reports during final submission checks.

## Author Contributions

Qiyun Hu: Conceptualization, Data curation, Formal analysis, Methodology, Visualization, and Writing – original draft. Ziqi Wang: Data curation and Formal analysis. Weina Wang: Formal analysis and Writing – review & editing. Hu Zhang: Conceptualization, Methodology, Supervision, and Writing – review & editing.

Final author approval status: `FINAL_AUTHOR_APPROVAL_REQUIRED`.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Conflict of Interest

The authors declare no conflicts of interest.

## Acknowledgments

[ACKNOWLEDGMENTS, IF ANY, TO BE PROVIDED. AI-ASSISTED WORKFLOW AND DRAFTING SHOULD BE DISCLOSED USING THE APPROVED TEXT IN `AI_use_disclosure_draft.md` AND THE JOURNAL'S FINAL PLACEMENT REQUIREMENT.]

## References

El Aidy, S., Bolsius, Y. G., Raven, F., & Havekes, R. (2020). A brief period of sleep deprivation leads to subtle changes in mouse gut microbiota. *Journal of Sleep Research, 29*(6), e12920. https://doi.org/10.1111/jsr.12920

Gao, T., Wang, Z., Dong, Y., Cao, J., Lin, R., Wang, X., Yu, Z., & Chen, Y. (2019). Role of melatonin in sleep deprivation-induced intestinal barrier dysfunction in mice. *Journal of Pineal Research, 67*(1), e12574. https://doi.org/10.1111/jpi.12574

Germain, P.-L., Lun, A., Garcia Meixide, C., Macnair, W., & Robinson, M. D. (2022). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research, 10*, 979. https://doi.org/10.12688/f1000research.73600.2

Karl, J. P., Whitney, C. C., Wilson, M. A., Fagnant, H. S., Radcliffe, P. N., Chakraborty, N., Campbell, R., Hoke, A., Gautam, A., Hammamieh, R., & Smith, T. J. (2023). Severe, short-term sleep restriction reduces gut microbiota community richness but does not alter intestinal permeability in healthy young men. *Scientific Reports, 13*, 213. https://doi.org/10.1038/s41598-023-27463-0

Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2021). Fast gene set enrichment analysis. *bioRxiv*, 060012. https://doi.org/10.1101/060012

Lee, J., Kang, J., Kim, Y., Lee, S., Oh, C.-M., & Kim, T. (2023). Integrated analysis of the microbiota-gut-brain axis in response to sleep deprivation and diet-induced obesity. *Frontiers in Endocrinology, 14*, 1117259. https://doi.org/10.3389/fendo.2023.1117259

Li, G., Gao, M., Zhang, S., Dai, T., Wang, F., Geng, J., Rao, J., Qin, X., Qian, J., Zuo, L., Zhou, M., Liu, L., & Zhou, H. (2024). Sleep deprivation impairs intestinal mucosal barrier by activating endoplasmic reticulum stress in goblet cells. *The American Journal of Pathology, 194*(1), 85–100. https://doi.org/10.1016/j.ajpath.2023.10.004

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database hallmark gene set collection. *Cell Systems, 1*(6), 417–425. https://doi.org/10.1016/j.cels.2015.12.004

Park, Y. S., Kim, S. H., Park, J. W., Kho, Y., Seok, P. R., Shin, J.-H., Choi, Y. J., Jun, J.-H., Jung, H. C., & Kim, E. K. (2020). Melatonin in the colon modulates intestinal microbiota in response to stress and sleep deprivation. *Intestinal Research, 18*(3), 325–336. https://doi.org/10.5217/ir.2019.00093

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods, 14*(4), 417–419. https://doi.org/10.1038/nmeth.4197

Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics, 26*(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

Rodríguez-Colman, M. J., Schewe, M., Meerlo, M., Stigter, E., Gerrits, J., Pras-Raves, M., Sacchetti, A., Hornsveld, M., Oost, K. C., Snippert, H. J., Verhoeven-Duif, N., Fodde, R., & Burgering, B. M. T. (2017). Interplay between metabolic identities in the intestinal crypt supports stem cell function. *Nature, 543*(7645), 424–427. https://doi.org/10.1038/nature21673

Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences. *F1000Research, 4*, 1521. https://doi.org/10.12688/f1000research.7563.2

Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck, W. M., III, Hao, Y., Stoeckius, M., Smibert, P., & Satija, R. (2019). Comprehensive integration of single-cell data. *Cell, 177*(7), 1888–1902.e21. https://doi.org/10.1016/j.cell.2019.05.031

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences of the United States of America, 102*(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102

Vaccaro, A., Kaplan Dor, Y., Nambara, K., Pollina, E. A., Lin, C., Greenberg, M. E., & Rogulja, D. (2020). Sleep loss can cause death through accumulation of reactive oxygen species in the gut. *Cell, 181*(6), 1307–1328.e15. https://doi.org/10.1016/j.cell.2020.04.049

Zhang, H.-Y., Shu, Y.-Q., Li, Y., Hu, Y.-L., Wu, Z.-H., Li, Z.-P., Deng, Y., Zheng, Z.-J., Zhang, X.-J., Gong, L.-F., Luo, Y., Wang, X.-Y., Li, H.-P., Liao, X.-P., Li, G., Ren, H., Qiu, W., & Sun, J. (2024). Metabolic disruption exacerbates intestinal damage during sleep deprivation by abolishing HIF1alpha-mediated repair. *Cell Reports, 43*(11), 114915. https://doi.org/10.1016/j.celrep.2024.114915

Zhang, M., Wu, X., Liu, D., Li, H., Li, X., Yang, W., Ye, J., Hou, L., Wang, S., Ning, N., Zhang, H., Tian, Y., Yu, L., Wu, K., Wang, L., Plikus, M. V., Lv, C., Wang, F., & Yu, Z. (2026). Sleep disturbance triggers aberrant activation of vagus circuitry and induces intestinal stem cell dysfunction. *Cell Stem Cell, 33*(2), 306–324.e8. https://doi.org/10.1016/j.stem.2026.01.002

## Figure Legends

### Figure 1. Study design, systematic dataset selection, and single-cell atlas

**(A)** Public-dataset selection and analysis sequence. Twenty-three unique candidate studies were assessed; three were eligible, one was potentially eligible, and 19 were excluded. GSE289089 was used for discovery, followed by independent validation in PRJNA1167170 and PRJNA862187. **(B)** UMAP of 27,671 GSE289089 singlets colored by 18 frozen clusters. **(C)** UMAP colored by manuscript cell-group mapping. **(D)** Canonical-marker DotPlot reconstructed from the frozen object. **(E)** Descriptive sample-level proportions of primary epithelial groups; no compositional inference was performed.

### Figure 2. Cell-type heterogeneity of oxidative-phosphorylation transcriptional responses

**(A)** Hallmark OXPHOS NES values for eight primary epithelial groups. Symbols denote FDR below 0.05. **(B)** NES matrix for Hallmark OXPHOS and five prespecified secondary mitochondrial pathways. **(C–F)** Running enrichment-score curves for intestinal stem cells, transit-amplifying progenitors, tuft cells, and serotonergic enterochromaffin cells.

### Figure 3. Independent bulk validation across two eligible colon datasets

**(A)** PCA of PRJNA1167170 proximal-colon samples. **(B)** PRJNA1167170 primary Hallmark OXPHOS enrichment. **(C)** PCA of the 12 PRJNA862187 colon samples. **(D)** PRJNA862187 primary average-effect Hallmark OXPHOS enrichment. **(E)** Three-dataset comparison of Hallmark OXPHOS and five prespecified secondary mitochondrial pathway NES values. Stars indicate FDR below 0.05; tissue-level datasets and intestinal-stem-cell pseudobulk differ in biological resolution.

### Figure 4. Cross-dataset gene-level concordance and direction

**(A)** GSE289089 intestinal-stem-cell versus PRJNA1167170 proximal-colon OXPHOS gene log-fold changes. Representative genes selected across mitochondrial modules for visualization are labelled. **(B)** Three-dataset directional summary for the 191 Hallmark genes measurable in all datasets, including the 82 genes negative in all three. Detailed pairwise correlations involving PRJNA862187 are reported in the text and Supporting Information.

### Figure 5. Robustness and sensitivity analyses

**(A)** Hallmark OXPHOS NES values under log-fold-change and signed-statistic rankings. **(B)** NES distributions across 20 equal-cell downsampling iterations for six cell groups. **(C)** Descriptive sample-level single-cell quality metrics; no inferential test is shown.

## Supporting Information legend

**Figure S1.** Single-cell quality-control distributions. **Figure S2.** Canonical-marker support for frozen annotations. **Figure S3.** Sample-level cell-group distributions. **Figure S4.** PRJNA1167170 Salmon diagnostics. **Figure S5.** PRJNA862187 Salmon diagnostics and sample metadata. **Figure S6.** PRJNA862187 diet-specific OXPHOS enrichment. **Figure S7.** PRJNA862187 exploratory interaction results. **Figure S8.** Mitochondrial pathway enrichment across primary and exploratory cell groups. **Figure S9.** Equal-cell downsampling under both ranking methods. **Figure S10.** Mixed regeneration-associated programs. **Figure S11.** Systematic public-dataset selection flow.

**Table S1.** Public-dataset candidate registry and exclusion reasons. **Table S2.** PRJNA862187 sample metadata and Salmon QC. **Table S3.** PRJNA862187 primary, diet-specific, secondary-pathway, and interaction results. **Table S4.** Three-dataset Hallmark OXPHOS gene-level concordance. **Table S5.** Three-dataset mitochondrial pathway summary.

---

**Draft status:** source-integrated submission draft. Repository deposition and final author approval remain mandatory before submission.
