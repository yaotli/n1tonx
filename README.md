1. Geographical distribution of Gs/GD viruses based on available HA sequences (Fig. Geo)
   - (a) 
     - ML by IQ-Tree
     - HA sequences (n=8487)
     - Curated from raw data with 3 outliers removed (`raw_seq.R`, `GsGD.R`)
   - (b)(c)
     - Clade classification based on the same tree as Fig. Geo(a)
     - Clade 2.3.2 (n=1802); clade 2.3.4-N1 (n=637); clade 2.3.4.4 (n=3214)

2. MCC tree of Gs/GD H5N1 lineages in China (Fig. GsGD-MCC)
   - MCMC run and annotated by BEAST (`HA_gsgd_sam3.xml`)
   - HA sequences (n=369)
   - Isolates in China were extracted from the alignment used in Fig. Geo(a), and the duplicated sequences were removed. Facilitated by an annotated tree labeled with (1) reference strains (2) groups containing viruses isolated in the same host or outbreak (`ha-gsgd.tre`), the alignment was subsampled by randomly selecting 25 isolate/group each year (`GsGD.R`, `_gsgd_subsample_.R`). Outliers were detected and removed before analyses by TempEst

3. Characterization of GsGD viruses in China based on available H5-HA sequences (Fig. Seq_info)
   - (a)(b)
     - HA sequences
     - Isolates in China were extracted from the alignment used in Fig. Geo(a), and clades were identified in sequence set with duplicated sequences (n=2480) and without duplicated sequences (n=1800)
   - (c)
     - States were determined by metadata appended with sequence in the database or past literatures (`_ecology_state.R`, `eco_states.tsv`)

4. Population dynamics of the three H5N1 lineages in China (Fig. pop_dynamic & a detailed figure)
   - (a)
     - MCMC run with Skygrid prior by BEAST (example files: `HA_c232_sam2.xml`...)
     - Sequences of the three clades were extracted from the alignment used in Fig. Geo(a) after constructing a IQ-tree (`H5_clade.R`). Duplicated sequences and sequences not isolated in China were then removed for each clade
     - Subsampling and identifying of outliers were applied as Fig. GsGD-MCC (`_clades_subsample_.R`). The resulting sequence set: clade 2.3.2 (n=122-123); clade 2.3.4 (n=201); clade 7 (n=55)
     - Growth rate by Skygrowth (https://github.com/mrc-ide/skygrowth)
   - (b)
     - Analyses as (a) (`N1_c232_sam2.xml`...)
     - Corresponding N1 sequences were first identified and processed as (a)(`N1_clade.R`). Clade 2.3.2 (n=135), clade 2.3.4 (n=88); clade 7 (n=35)

5. Population dynamics of the H6N1 lineage in China (Fig. Pop_h6)
     - Analyses as Fig. pop_dynamic(a) (`H6_sam2.xml`, `H6_scA2.xml`)
     - H6 alignment preparation was similar to H5/N1 (`nonH5.R`). Observed (n=182); sampling scheme (n=184)

5. Coevolution of H5 and N1 genes of GsGD viruses in China (Fig. HA-NA)
   - Left panel 
     - H5 tree (n=8487) is identical to the one used in Fig. Geo(a), whereas N1 tree was built by H5's corresponding N1 (n=5232)
   - Right panel
     - Generated along with `Fig. pop_dynamic(b)` analyses

6. Tests of constant population size based on nucleotide polymorphism (Fig. Pop_stat)
   - HA sequences
   - The sequence sets used for each clade is identical to the one prepared in `Fig. pop_dynamic(a)` before subsampling - clade 2.3.2 (n=276), clade 2.3.4 (n=1085), clade 7 (n=74)
   - Methods: Tajima's D & pi (Nei and Kumar, 2000); Fu and Li's D & F (Fu and Li 1993, Simonsen et al., 1995)

7. Table 
   - Clade 2.3.4.4 (internal branch)
     - Parsed in `Fig. pop_dynamic(a)` MCC trees
   - Clade 2.3.4.4 (root)
     - The alignment used here was prepared as clade 2.3.4 in `Fig. pop_dynamic(a)` but subsampling was only conducted with clade 2.3.4.4 (n=118)
   - Clade 7
     - Parsed in `Fig. pop_dynamic(a)` MCC trees

8. Gene flow between ecological systems (Fig. gene_flow)
   - MCMC run and annotated by BEAST (`HA_c232_eco.xml`...)
   - HA sequences
   - Sequences for each clade were identical to sequences subsampled in `Fig. pop_dynamic`. Wild or domestic states were designated with `eco_states.tsv`. Clade 2.3.2 (n=123), clade 2.3.4 (n=201); clade 7 (n=55). 
   - Trunk proportions by PACT (https://github.com/trvrb/PACT)

9. Transition parameters between the two ecological states (Fig. Eco_param)
   - Parameters parsed from Bayesian discrete diffusion model (MCC trees in `Fig. gene_flow`).
   - Another sets of alignments were prepared with an increased in maximum sampled sequences each year (25 -> 30). Clade 2.3.2 (n=122); clade 2.3.4 (n=223); clade (n=55). 

10. Epidemiology of Gs/GD viruses and poultry production in China (Fig. Epi)
    - Data source for number of outbreaks (http://empres-i.fao.org/)
    - Number of cases (https://www.who.int/influenza/human_animal_interface/H5N1_cumulative_table_archives/en/)
    - Poultry production in China (http://www.fao.org/faostat/en/#data)
	  
