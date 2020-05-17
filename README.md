1. Sup_fig 1
   - (a) 
     - ML by IQ-Tree
     - HA sequecnes (n=8487)
     - Curated from raw data with 3 outliers removed (`raw_seq.R`, `GsGD.R`)
   - (b)(c)
     - Clade classification based on the same tree as Sup_fig 1a
     - Clade 2.3.2 (n=1802); clade 2.3.4-N1 (n=637); clade 2.3.4.4 (n=3214)

2. Fig 1
   - MCMC run and annotated by BEAST (`HA_gsgd_sam3.xml`)
   - HA sequences (n=369)
   - Isolates in China were extracted from the alignment used in Sup_fig 1(a), and the duplicated sequences were removed. Facilitated by an annotated tree labeled with (1) reference strains (2) groups containing viruses isolated in the same host or outbreak (`ha-gsgd.tre`), the alignment was subsampled by randomly selecting 25 isolate/group each year (`GsGD.R`, `_gsgd_subsample_.R`). Outliers were detected and removed before analyses by TempEst

3. Sup_fig 2
   - (a)(b)
     - HA sequences
     - Isolates in China were extracted from the alignment used in Sup_fig 1(a), and clades were identified in sequence set with duplicated sequences (n=2480) and without duplicated sequences (n=1800)
   - (c)
     - States were determined by metadata appended with sequence in the database or past literatures (`_ecology_state.R`, `eco_states.tsv`)

4. Figure 2 and Sup_fig 3
   - (a)
     - MCMC run with Skygrid prior by BEAST (example files: `HA_c232_sam2.xml`...)
     - Sequences of the three clades were extracted from the alignment used in Sup_fig 1 after constructing a IQ-tree (`H5_clade.R`). Duplicated sequences and sequences not isolated in China were then removed for each clade
     - Subsampling and identifying of outliers were applied as Fig 1 (`_clades_subsample_.R`). The resulting sequence set: clade 2.3.2 (n=122-123); clade 2.3.4 (n=201); clade 7 (n=55)
     - Growth rate by Skygrowth (https://github.com/mrc-ide/skygrowth)
   - (b)
     - Analyses as (a) (`N1_c232_sam2.xml`...)
     - Corresponding N1 sequences were first identified and processed as (a)(`N1_clade.R`). Clade 2.3.2 (n=135), clade 2.3.4 (n=88); clade 7 (n=35)
   - (c)
     - Analyses as (a) (`H6_sam2.xml`, `H6_scA2.xml`)
     - H6 alignment preparation was similar as H5/N1 (`nonH5.R`). Observed (n=182); sampling scheme (n=184)

5. Sup_fig 4
   - Left panel 
     - H5 tree (n=8487) is identical to the one used in Sup_fig 1a, whereas N1 tree was built by H5's corresponding N1 (n=5232)
   - Right panel
     - Generated from `Figure 2b` analyses

6. Figure 3
   - HA sequences
   - The sequence sets used for each clade is identical to the one prepared in `Figure 2a` before subsampling - clade 2.3.2 (n=276), clade 2.3.4 (n=1085), clade 7 (n=74)
   - Methods: Tajima's D & pi (Nei and Kumar, 2000); Fu and Li's D & F (Fu and Li 1993, Simonsen et al., 1995)

7. Table 
   - Clade 2.3.4.4 (internal branch)
     - Parsed in `Figure 2a` MCC 
   - Clade 2.3.4.4 (internal branch)
     - The alignment was prepared as clade 2.3.4 in `Figure 2a` but subsampling was only conducted with clade 2.3.4.4 (n=118)
   - Clade 7
     - Parsed in `Figure 2a` MCC

8. Figure 4 and Supp_fig 5
   - MCMC run and annotated by BEAST (`HA_c232_eco.xml`...)
   - HA sequences
   - Sequences for each clade were identical to sequences subsampled in `Figure 2`. Wild or domestic states were designated with `eco_states.tsv`. Clade 2.3.2 (n=123), clade 2.3.4 (n=201); clade 7 (n=55). 
   - Trunk proportions by PACT (https://github.com/trvrb/PACT)

9. Figure 5
   - Parameters parsed from Bayesian discrete diffusion model (MCC trees in `Figure 4` and `Supp_fig 5`).
   - Another sets of alignments were prepared with an increased in maximum sampled sequences each year (25 -> 30). Clade 2.3.2 (n=122); clade 2.3.4 (n=223); clade (n=55). 

10. Supp_fig 6
    - Data source for number of outbreaks (http://empres-i.fao.org/)
    - Number of cases (https://www.who.int/influenza/human_animal_interface/H5N1_cumulative_table_archives/en/)
    - Poultry production in China (http://www.fao.org/faostat/en/#data)


	  
