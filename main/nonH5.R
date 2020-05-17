###### Raw extraction of H9N2/H6Nx ######

# This step is nearly identical to the one that applied to pH5Nx and subclade (& their N1). 

source( "../functions/f_strain.R" )

# 1 copy the files from `raw` to `processed` ------

# # cml1
# cp -p N1TONX/raw/pH9N2/align/rm_gap_amb/S4_pH9-rsa-rga-rga_trim.fasta N1TONX/processed/pH9N2/raw
# for f in $(ls *fasta); do mv "$f" "${f/-*.fasta/-slc.fasta}"; done
# cp -p N1TONX/raw/pH9N2/isolate_ac.2.tsv N1TONX/processed/pH9N2/raw
# cp -p N1TONX/raw/pH9N2/isolate_meta.tsv N1TONX/processed/pH9N2/raw

# # cml2
# cp -p N1TONX/raw/pH9N2/align/rm_gap_amb/S4_pH9-rsa-rga-rga_trim.fasta N1TONX/processed/pH9N2/raw
# for f in $(ls *fasta); do mv "$f" "${f/-*.fasta/-slc.fasta}"; done
# cp -p N1TONX/raw/pH6Nx/isolate_ac.tsv N1TONX/processed/pH6Nx/raw
# cp -p N1TONX/raw/pH6Nx/isolate_meta.tsv N1TONX/processed/pH6Nx/raw


# 2 update the sequence "name" ------

# h9
# slc_fasta   = "N1TONX/processed/pH9N2/raw/"
# slc_fasta_i = list.files( slc_fasta, full.names = TRUE, pattern = "-slc.fasta" )
# ac_tab_xi   = read.table( paste0( slc_fasta, "isolate_ac.2.tsv"),
#                           header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )
# 
# for( i in 1: length(slc_fasta_i) )
# {
#   id0  = fastaEx( slc_fasta_i[i] )$id
#   seq0 = fastaEx( slc_fasta_i[i] )$seq
#   seg  = as.numeric( str_match( slc_fasta_i[i], "S([0-9])_" )[,2] )
# 
#   id_out = map( id0,
#                 function(x)
#                 {
#                   m_i = match( x, ac_tab_xi[, seg] )
# 
#                   y =
#                   ifelse( !is.na(m_i),
#                           paste0( rownames(ac_tab_xi)[m_i], "-", x ),
#                           paste0( "UNID", "-", x )  )
#                   return(y)
#                 })
# 
#   write.fasta( seq0, id_out, gsub( ".fasta", "-x1.fasta", slc_fasta_i[i]) )
#   print(i)
# }
# 
# ls_slc_x1 = list.files( "N1TONX/processed/pH9N2/raw/", pattern = "slc-x1.fasta$", full.names = TRUE )
# ls_tsv    = list.files( "N1TONX/raw/pH9N2/", pattern = "S[0-9].*tsv$", full.names = TRUE )
# ls_tsv    = ls_tsv[1]
# 
# leaf( infile     = ls_slc_x1,
#       tsvfile    = ls_tsv,
#       includeOld = TRUE,
#       isTree     = FALSE,
#       include    = c( "name", "sero", "country", "year" ))

# h6
# slc_fasta   = "N1TONX/processed/pH6Nx/raw/"
# slc_fasta_i = list.files( slc_fasta, full.names = TRUE, pattern = "-slc.fasta" )
# ac_tab_xi   = read.table( paste0( slc_fasta, "isolate_ac.tsv"),
#                           header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )
# 
# for( i in 1: length(slc_fasta_i) )
# {
#   id0  = fastaEx( slc_fasta_i[i] )$id
#   seq0 = fastaEx( slc_fasta_i[i] )$seq
#   seg  = as.numeric( str_match( slc_fasta_i[i], "S([0-9])_" )[,2] )
#   
#   id_out = map( id0,
#                 function(x)
#                 {
#                   m_i = match( x, ac_tab_xi[, seg] )
#                   
#                   y =
#                     ifelse( !is.na(m_i),
#                             paste0( rownames(ac_tab_xi)[m_i], "-", x ),
#                             paste0( "UNID", "-", x )  )
#                   return(y)
#                 })
#   
#   write.fasta( seq0, id_out, gsub( ".fasta", "-x1.fasta", slc_fasta_i[i]) )
#   print(i)
# }
# 
# ls_slc_x1 = list.files( "N1TONX/processed/pH6Nx//raw/", pattern = "slc-x1.fasta$", full.names = TRUE )
# ls_tsv    = list.files( "N1TONX/raw/pH6Nx/", pattern = "S[0-9].*tsv$", full.names = TRUE )
# ls_tsv    = ls_tsv[1]
# 
# leaf( infile     = ls_slc_x1,
#       tsvfile    = ls_tsv,
#       includeOld = TRUE,
#       isTree     = FALSE,
#       include    = c( "name", "sero", "country", "year" ))
  

# fasttree

# 3 manual select outliers in each gene and remove ------
# no apparent outliers are spotted at these two tree


# 4 extract CN-HK sequences ------

source( "../functions/function.R" )

# path_tsv    = c( "N1TONX/raw/pH9N2/S4_pH9.tsv", "N1TONX/raw/pH6Nx/S4_pH6.tsv" )
# path_slc_x1 = c( "N1TONX/processed/pH9N2/slc-x1/S4_pH9-slc-x1.fasta", 
#                  "N1TONX/processed/pH6Nx/slc-x1/S4_pH6-slc-x1.fasta" )
# 
# df_tsv = do.call( rbind,
#                   lapply( as.list( path_tsv ),
#                           function(x) read.table( x, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) ) )
# 
# for( i in 1: length(path_slc_x1) )
# {
#   fas_seq = fastaEx( path_slc_x1[i] )$seq
#   fas_id  = fastaEx( path_slc_x1[i] )$id
# 
#   x =
#   map( fas_id,
#        function(x)
#        {
#          ac0 = str_split( x, "__")[[1]][1]
# 
#          ac = str_extract( ac0, "[A-Z0-9]+$" )
#          m  = match( ac, df_tsv$ac )
# 
#          if( is.na(m) | length(m) > 1 ){ stop("mismatch") }
# 
#          y = grepl( "China|Hong_Kong", df_tsv$country[m], ignore.case = TRUE )
#          return(y)
#        })
# 
#   z = unlist(x)
#   write.fasta( sequences = fas_seq[z], names = fas_id[z],
#                file.out  = gsub( "x1.fasta", "CN.fasta", path_slc_x1[i]  ) )
# }


# 5 remove duplicate sequences ------

source( "../functions/f_rmdup.R" )

# fast_rmdup_v3( "N1TONX/processed/pH9N2//slc-CN/S4_pH9-slc-CN.fasta",
#                "N1TONX/raw/pH9N2/S4_pH9.tsv" )
# 
# fast_rmdup_v3( "N1TONX/processed/pH6Nx//slc-CN/S4_pH6-slc-CN.fasta",
#                "N1TONX/raw/pH6Nx/S4_pH6.tsv" )


# 6 subsample Y280 lineage of H9N2 and group2 (ST2853-like) of H6 ------

source( "_clades_subsample_.R" )

# check 

# tre_df_com = data.frame()
# for( i in 1:2 )
# {
#   grouped_trees = c( "N1TONX/processed/pH6Nx/slc-CN/rmdup/iqtree/S4_pH6-slc-CN_rmDup_col4.tre",
#                      "N1TONX/processed/pH9N2//slc-CN/rmdup/iqtree/S4_pH9-slc-CN-Y280_col2.tre" )
#                      
#   crude_nex = read.csv( grouped_trees[i], stringsAsFactors = FALSE )
# 
#   taxa.s = grep( x = crude_nex[,1], pattern = "taxlabels" ) + 1
#   ntax   = as.numeric( str_match( grep( "ntax", crude_nex[,1], value = TRUE), "(ntax=)([0-9]+)")[,3] )
#   taxa.e = taxa.s + ntax - 1
# 
#   tre_name  = gsub( "\t|\\[.*\\]", "", crude_nex[, 1][taxa.s: taxa.e] )
#   tre_col   = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "color=#([0-9a-z]+)" )[,2]
#   tre_group = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "name=([0-9]+)" )[,2]
# 
#   tre_df = data.frame( name = tre_name, col = tre_col, group = tre_group, stringsAsFactors = FALSE )
# 
#   print(i)
#   tre_df_com = rbind( tre_df_com, tre_df )
# }
# 
# login = c()
# for( j in 1: length( unique( tre_df_com$group ) ) )
# {
#   gp = unique( tre_df_com$group )[j]
#   if( is.na(gp) ){ next() }
# 
#   names_i = tre_df_com$name[ which( tre_df_com$group == gp ) ]
#   print( names_i )
# 
#   keyin = readline( prompt = "1 for sth. wrong else enter" )
# 
#   if( keyin == "k" ){ stop() }
# 
#   login = c( login, keyin )
# 
#   print( paste0( j, " / ", length( unique( tre_df_com$group ) )) )
# }

# pre_sumsample( "N1TONX/processed/pH6Nx/slc-CN/rmdup/iqtree/S4_pH6-slc-CN_rmDup_col4.tre",
#                yr_range     = 1,
#                batchname    = "H6Nx", 
#                col_toremove = "bcbd22" )
# 
# pre_sumsample( "N1TONX/processed/pH9N2//slc-CN/rmdup/iqtree/S4_pH9-slc-CN-Y280_col2.tre" ,
#                yr_range  = 3,
#                batchname = "H9N2")
# pre_sumsample( "N1TONX/processed/pH9N2//slc-CN/rmdup/iqtree/S4_pH9-slc-CN-Y280-c1_col2.tre" ,
#                yr_range  = 3,
#                batchname = "H9N2-c1")
# pre_sumsample( "N1TONX/processed/pH9N2//slc-CN/rmdup/iqtree/S4_pH9-slc-CN-Y280-c2_col2.tre" ,
#                yr_range  = 3,
#                batchname = "H9N2-c2")


# subsample 

# input_sam_table = c( "N1TONX/processed/pH6Nx/CN-x1/H6Nx_sam_table.tsv",
#                      "N1TONX/processed/pH9N2/CN-x1/H9N2_sam_table.tsv",
#                      "N1TONX/processed/pH9N2/CN-x1/H9N2-c1_sam_table.tsv",
#                      "N1TONX/processed/pH9N2/CN-x1/H9N2-c2_sam_table.tsv" )
# 
# input_fas       = c( "N1TONX/processed/pH6Nx/slc-CN/rmdup/S4_pH6-slc-CN_rmDup.fasta",
#                      "N1TONX/processed/pH9N2/slc-CN/rmdup/S4_pH9-slc-CN_rmDup.fasta" )
# 
# fas_id = unlist( sapply( as.list( input_fas ),
#                           function(x)
#                           {
#                             y = fastaEx(x)$id
#                           } ) )
# 
# fas_seq = lapply( as.list( input_fas ),
#                          function(x)
#                          {
#                            y = fastaEx(x)$seq
#                          } )
# fas_seq = do.call( c, fas_seq )
# 
# nseed = c( 123, 234, 345, 456, 567 )
# 
# for( i in 1: length(input_sam_table) )
# {
#   tab = read.table( input_sam_table[i], sep = "\t", stringsAsFactors = FALSE, header = TRUE )
# 
#   for( j in 1:5 )
#   {
#     ran = subsample_gsgd( xseed = nseed[j], path_sam_table = input_sam_table[i] )
# 
#     m = match( tab$name[ ran ], fas_id )
#     if( TRUE %in% is.na(m) ){ stop( "mismatch" ) }
# 
#     if( i == 1 ){ filename = paste0( "N1TONX/processed/pH6Nx/CN-x1/S4_pH6_sam", j, ".fasta" ) }
#     else if( i == 2 )
#     {
#       filename = paste0( "N1TONX/processed/pH9N2/CN-x1/S4_pH9_sam", j, ".fasta" )
#     }
#     else if( i == 3 )
#     {
#       filename = paste0( "N1TONX/processed/pH9N2/CN-x1/S4_pH9-c1_sam", j, ".fasta" )
#     }
#     else
#     {
#       filename = paste0( "N1TONX/processed/pH9N2/CN-x1/S4_pH9-c2_sam", j, ".fasta" )  
#     }
# 
#     write.fasta( fas_seq[m], fas_id[m], filename )
#   }
#  print(i)
# }

# 7 remove outliers ------

# all-H9N2

# path_ls_fas  = list.files( "N1TONX/processed/pH9N2/CN-x1/", pattern = ".fasta", full.names = TRUE )
# path_outlier = "N1TONX/processed/pH9N2/CN-x2/outliers.txt"
# 
# toberemoved = read.table( path_outlier, sep = "\t", stringsAsFactors = FALSE, header = FALSE )
# toberemoved = unique( toberemoved$V1 )
# 
# for( f in 1:length( path_ls_fas ) )
# {
#   seq_id  = fastaEx( path_ls_fas[f] )$id
#   seq_fas = fastaEx( path_ls_fas[f] )$seq
#   
#   m = match( toberemoved, seq_id )
#   m = m[ !is.na(m) ]
#   
#   write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x2.fasta", path_ls_fas[f] ) ) 
# }

# H9N2-c1/2

# path_ls_fas  = c( list.files( "N1TONX/processed/pH9N2/CN-x1/", pattern = "pH9-c1", full.names = TRUE ),
#                   list.files( "N1TONX/processed/pH9N2/CN-x1/", pattern = "pH9-c2", full.names = TRUE ) )
# path_outlier = c( "N1TONX/processed/pH9N2/CN-c1/outliers.txt", 
#                   "N1TONX/processed/pH9N2/CN-c2/outliers.txt" )
# 
# toberemoved = do.call(rbind, lapply( as.list( path_outlier ),
#                                      function(x)
#                                      {
#                                        read.table( x, sep = "\t", stringsAsFactors = FALSE, header = FALSE )
#                                      }) )
# toberemoved = unique( toberemoved$V1 )
# 
# for( f in 1:length( path_ls_fas ) )
# {
#   seq_id  = fastaEx( path_ls_fas[f] )$id
#   seq_fas = fastaEx( path_ls_fas[f] )$seq
# 
#   m = match( toberemoved, seq_id )
#   m = m[ !is.na(m) ]
# 
#   write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x2.fasta", path_ls_fas[f] ) )
# }


# all-H6

# path_ls_fas  = list.files( "N1TONX/processed/pH6Nx/CN-x1/", pattern = ".fasta", full.names = TRUE )
# path_outlier = "N1TONX/processed/pH6Nx/CN-x2/outliers.txt"
# 
# toberemoved = read.table( path_outlier, sep = "\t", stringsAsFactors = FALSE, header = FALSE )
# toberemoved = unique( toberemoved$V1 )
# 
# for( f in 1:length( path_ls_fas ) )
# {
#   seq_id  = fastaEx( path_ls_fas[f] )$id
#   seq_fas = fastaEx( path_ls_fas[f] )$seq
# 
#   m = match( toberemoved, seq_id )
#   m = m[ !is.na(m) ]
# 
#   write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x2.fasta", path_ls_fas[f] ) )
# }


###### Sampling schemes for H6 population ######

source( "_clades_subsample_.R" )

# input_sam_table = "N1TONX/processed/pH6Nx/CN-x1/H6Nx_sam_table.tsv"
# input_fas       = "N1TONX/processed/pH6Nx/slc-CN/rmdup/S4_pH6-slc-CN_rmDup.fasta" 
# 
# fas_id  = fastaEx( input_fas )$id
# fas_seq = fastaEx( input_fas )$seq
# 
# nseed = c( 123, 234, 345, 456, 567 )
# 
# for( i in 1: length(input_sam_table) )
# {
#   tab = read.table( input_sam_table[i], sep = "\t", stringsAsFactors = FALSE, header = TRUE )
# 
#   for( j in 1:5 )
#   {
#     ran = subsample_gsgd_v2( xseed = nseed[j], path_sam_table = input_sam_table[i], 
#                              max_n       = 30, 
#                              yr_resample = c(2007,2008),
#                              yr_prop     = c(0,0.5) )
# 
#     m = match( tab$name[ ran ], fas_id )
#     
#     if( TRUE %in% is.na(m) ){ stop( "mismatch" ) }
# 
#     filename = paste0( "N1TONX/processed/pH6Nx/CN-x1/S4_pH6_scA", j, ".fasta" )
#     
#     write.fasta( fas_seq[m], fas_id[m], filename )
#   }
# }

# outlier

path_ls_fas  = list.files( "N1TONX/processed/pH6Nx/SC-a1/", pattern = ".fasta", full.names = TRUE )
path_outlier = "N1TONX/processed/pH6Nx/SC-a2/outlier.txt"

toberemoved = read.table( path_outlier, sep = "\t", stringsAsFactors = FALSE, header = FALSE )
toberemoved = unique( toberemoved$V1 )

for( f in 1:length( path_ls_fas ) )
{
  seq_id  = fastaEx( path_ls_fas[f] )$id
  seq_fas = fastaEx( path_ls_fas[f] )$seq

  m = match( toberemoved, seq_id )
  m = m[ !is.na(m) ]

  write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x2.fasta", path_ls_fas[f] ) )
}


