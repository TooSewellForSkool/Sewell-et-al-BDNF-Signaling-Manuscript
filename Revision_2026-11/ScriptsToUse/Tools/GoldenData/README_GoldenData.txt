GOLDEN DATA FOR THE STATISTICS CHECK
=====================================
PublishedStats_TableS3_long.csv
  What:    the 220 p-values of Supplementary Table 3 (Fig 2C and Supplementary Fig 4B), one row per
           marker x time point x treatment, written in long form so Compare_GoldenFiles.R can match rows.
  Source:  Supplementary_Tables.xlsx, sheet 'Table S3' (file MD5 fccc9045db04a98304ba702d92951feb), sent by JMS on 2026-10-25.
  Made by: reading the sheet cell by cell; the only change is that the column headings used in the
           paper (pERK, pAkt, pCREB, ...) were replaced by the pipeline marker names (pERK_171Yb, ...).
           Values are copied as they are in the sheet (the sheet shows 7-10 significant digits).
  Why a looser tolerance: the sheet was assembled from printed tables, so values are rounded;
           the check uses 1e-7 (the largest difference found when the p-values were recomputed
           independently in Python from the published RAW table was 4.5e-8).
  Groups:  NoTreat (time point 0, p = 1 because it is compared with itself), ContDepr, BDNF.

PublishedStats_TableS4A_long.csv, ..._TableS4B_long.csv, ..._TableS4C_long.csv
  What:    the p-values of Supplementary Table 4 panels A (Fig 3C: clusters 1,2,4,7,13; 150 values),
           B (Fig 3D: clusters 1,3,15,17; 96 values) and C (Fig 3E: clusters 9,17,18; 90 values), long form:
           one row per marker x time point x cluster (Group = cluster number).
  Source:  Supplementary_Tables.xlsx, sheet 'Table S4' (file MD5 fccc9045db04a98304ba702d92951feb); same conversion as Table S3.
           The treatment compared with NoTreat is not named in the sheet; it is taken from the published folder
           names (Fig 3C and 3E: BDNF; Fig 3D: CompMed).
  Note:    at time point 0 every cluster has p = 1: the fold-change table is normalised so every cluster's 0 hr mean is 1,
           and the published code compared each cluster's 0 hr values with the 0 hr values pooled over the plotted clusters.
           That pattern shows the published tests were run on the FOLD CHANGE table (not on percent above threshold).
