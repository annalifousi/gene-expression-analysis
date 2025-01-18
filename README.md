Low-Dose Radiation Therapy and CAR T-Cell Efficacy in CD19+ Lymphoma

This project investigates the impact of low-dose radiation therapy (RT) on CAR T-cell efficacy in a mouse model of CD19+ lymphoma, focusing on gene expression changes in irradiated (IR) vs. non-irradiated (Non-IR) tumors and lymph nodes at 24 hours and 7 days. The goal is to uncover immune pathways and tumor responses that can inform personalized treatments for hematologic malignancies.


Data Avaible at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281695


Objectives
Gene Expression Analysis: Investigate changes in gene expression between IR and Non-IR tumors and lymph nodes at two time points.
Human Precision Medicine: Identify cross-species gene homologs to guide personalized therapies, optimize dosing, and improve outcomes in hematologic cancers.


Data Analysis Overview

1. PCA Plot
Visualizes sample grouping based on IR vs. Non-IR conditions and time points (24h vs. 7d).

2. Volcano Plots
Identifies significantly up- or downregulated genes between conditions.

3. Barplot of Enriched Pathways
Shows enriched biological pathways based on differentially expressed genes (DEGs), with focus on immune-related pathways.

4. Heatmap of Temporal Dynamics
Displays expression patterns of key DEGs over time and conditions.

==DESeq Analysis==

Purpose: Perform differential expression analysis between IR vs. Non-IR conditions using DESeq2.
Method: DESeq2 is used to identify DEGs, and results are visualized in volcano plots and heatmaps. This analysis focuses on determining the genes whose expression is significantly altered by radiation.

==Gene Set Enrichment Analysis (GSEA)==

Purpose: Identify enriched biological pathways using ranked DEGs.
Method: GSEA tests pathway enrichment by comparing the ranking of DEGs between conditions.


==EdgeR Analysis==

Identifies genes with significant expression changes due to radiation therapy. These DEGs are mapped to human orthologs to assess potential therapeutic relevance.
