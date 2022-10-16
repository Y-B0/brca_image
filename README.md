# brca_image
Correlate Histopathological Image and Transcriptomics Data to Identify the Molecular Basis of Nuclear Morphology with Gene Co-expression Network Approach


final_script.R: R script of project

TCGA.BRCA.htseq_counts.tsv.gz: BRCA read counts from GDC

gencode.v22.annotation.gene: Gene information file

BRCA.im.select.txt: morphology features value

genecount.txt: source read counts input file of TSUNAMI

raw_genecount.RData: processed read counts file of RData

protein_coding.hg38.position:mRNA position information

coexp_modules.csv: gene sets of modules from TSUNAMI

coexp_eigengene_matrix.csv: eigengene of modules from TSUNAMI

datamatrix_2.csv: read counts filtered by TSUNAMI

TCGA-BRCA.survival.tsv: survival data of BRCA from GDC

clinical.tsv: clinical information of BRCA from GDC

basic_data.RData: RData of all basically data

expr.txt: read counts data used to calculate immune score

Estimate_gene.gct: filtered read counts data from expr.txt

Estimate_score.gct: result of immune score

c1.all.v7.4.entrez.gmt: Gene sets corresponding to each human chromosome and each cytogenetic band

cytoband.txt: cytoband enrichment result

Table_GO_result.csv: GO term enrichment result 

GO_enrichment.RData: GO term enrichment result of RData file

enrich_model.RData: summary GO term enrichment result of all module 

cytoband_module.RData: summary cytoband enrichment result of all module

module*.tsv: tf enrichment result of modules

module*_regulate_network_allgene.txt: TF-gene network information of modules

module*_regulate_network_type_allgene.txt: node information of TF-gene network

module*_image_gene_heatmap.pdf: correlation heatmap of morphology feature and genes

Module*_interactions_short.tsv: ppi network file from STRING database

module*_cytohubba.csv: hubgene lists from ppi network

riskscore.RData: risk score result of samples

alter_predict_final.RData: predict model process file

multi_cox_result.txt: multivariate cox regression result of risk score 

ref.bib: bibtex of reference

module*_ppi.pdf: ppi network figure

module*_regulate_network.pdf: regulate network figure

univariate_cox_result.txt: univariate cox regression result of morphology features

Module*_expression.pdf: gene expression distribution of high- and low-risk score group 

feature_heatmap: heatmap of morphology features in high- and low-risk score group

feature_wilcox.pdf: morphology features wilcox test of high- and low-risk score group 

auc.pdf: AUC curve figure

survival.pdf: survival curve of high- and low-risk score group

riskscore.pdf: samples risk score distribution of high- and low-risk score group

feature_cor_eigengene.pdf: correlation heatmap of morphology features and eigengene

cytoband.pdf: cytobands enrichment result figure

go.pdf: go term enrichment result figure


![image](https://user-images.githubusercontent.com/52906974/196042952-55280eed-9086-4f06-84f5-c25657d26913.png)
