###### Raw extraction of three clades ######

# Three clades are firstly identified from `slc-x2-gsgd_colclade1.tre`, and then viruses 
# isolated in China (and Hong Kong) are extracted. As pH5Nx, the alignments undergo 
# rmdup(v3). With the rm_dup iqtrees, subsampling process is done following the procedure
# of pH5Nx. 

require( purrr )
require( stringr )

# 1 extract each clade ------

source("../functions/function.R")

# sed -i "" "s/,//g" S4_pH5-slc-x2-gsgd_colclade1.tre

# slc_x2_gsgd_colclade = "N1TONX/processed/pH5Nx/slc-x2/iqtree/S4_pH5-slc-x2-gsgd_colclade1.tre"
# raw_fas              = "N1TONX/processed/pH5Nx/slc-x2/S4_pH5-slc-x2.fasta"
# 
# tag_table = tagExtra( slc_x2_gsgd_colclade )
# 
# assigned_name = c( "c234", "c2344", "c234N1", "c232", "c7" )
# col_code      = c( "ff6666|ffcc66", "ffcc66", "ff6666", "ccff66", "66ccff" )
# 
# fas_id0  = fastaEx( raw_fas )$id
# fas_seq0 = fastaEx( raw_fas )$seq
# 
# for(i in 1:5)
# {
#   idx = grep( col_code[i], tag_table$tag )
#   
#   m = match( tag_table$id[idx], fas_id0 )
#   if( TRUE %in% is.na(m) ){ stop( "mismatch" ) }
#   
#   outname = gsub( "slc-x2.fasta", paste0( assigned_name[i], ".fasta" ), raw_fas )
#   
#   print( length(m) )
#   write.fasta( fas_seq0[m], fas_id0[m], outname )
# }

# 2 CN-HK ------

source( "../functions/function.R" )

# path_tsv     = "N1TONX/raw/pH5Nx/"
# path_sub_raw = "N1TONX/processed/H5_subclade/raw/"
# 
# input  = list.files( path_tsv, full.names = TRUE, pattern = "^S.*tsv$" )
# df_tsv = do.call( rbind,
#                   lapply( as.list( input ),
#                           function(x) read.table( x, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) ) )
# 
# list_fasta = list.files( path_sub_raw, pattern = ".fasta", full.names = TRUE )
# 
# for( i in 1: length(list_fasta) )
# {
#   fas_seq = fastaEx( list_fasta[i] )$seq
#   fas_id  = fastaEx( list_fasta[i] )$id
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
#                file.out  = gsub( ".fasta", "-CN.fasta", list_fasta[i]  ) )
# }

# 3 rmdup ------

source( "../functions/f_rmdup.R" )

# sub_CN   = list.files( "N1TONX/processed/H5_subclade/CN-x1/", pattern = "^S.*.fasta", full.names = TRUE)
# path_tsv = "N1TONX/raw/pH5Nx/S4_pH5.tsv"
# 
# for( i in 1:length(sub_CN) )
# {
#   fast_rmdup_v3( fas_file   = sub_CN[i],
#                  table_file = path_tsv )
# }


# 4 grouping ------

# manually group as gsgd dataset 

## check 
# tre_df_com = data.frame()
# for( i in 1:3 )
# {
#   grouped_trees = c( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c232-CN_rmDup_col2.tre",
#                      "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c234-CN_rmDup_col2.tre",
#                      "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c7-CN_rmDup_col2.tre" )
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

# 5 subsample ------

cat( "_clades_subsample_.R" )

# remove two outliers found in TempEst 

# toberemoved = read.table( "N1TONX/processed/H5_subclade/CN-x3/outlier.tsv", 
#                           sep = "\t", stringsAsFactors = FALSE, header = FALSE )
#   
# inputfas    = list.files( "N1TONX/processed/H5_subclade/CN-x2/", 
#                           pattern = ".fasta", full.names = TRUE )
# 
# for( f in 1:length(inputfas) )
# {
#   clade = str_match( inputfas[f], "c([0-9N]+)_sam[0-9]" )[,2]
#   
#   rr = which( toberemoved$V1 == clade )
#   
#   if( length(rr) == 0 )
#   {
#     cpname = gsub( ".fasta", "-x3.fasta", inputfas[f] )
#     cmd    = paste0( "cp ", inputfas[f], " ", cpname )
#     system( cmd )
#     
#   }else
#   {
#     seq_id  = fastaEx( inputfas[f] )$id
#     seq_fas = fastaEx( inputfas[f] )$seq  
#    
#     m = match( toberemoved$V2[rr], seq_id ) 
#     m = m[ !is.na(m) ]
#     
#     write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x3.fasta", inputfas[f] ) )
#   }
#   
# }

###### 8 BEAST xml ######

cat("_BEAST_xml_file_.R")

###### 11 prepare larger data set ######

source( "_clades_subsample_.R" )

# subsample 

# input_sam_table = list.files( "N1TONX/processed/H5_subclade/CN-x2",
#                               pattern = "sam_table", full.names = TRUE )[c(1,2,5)]
# input_fas       = list.files( "N1TONX/processed/H5_subclade/CN-x1",
#                               pattern = "^S[0-9]_.*.fasta$", full.names = TRUE )
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
# nseed = c( 123, 234, 345 )
# 
# for( i in 1: length(input_sam_table) )
# {
#   tab = read.table( input_sam_table[i], sep = "\t", stringsAsFactors = FALSE, header = TRUE )
# 
#   for( j in 1:3 )
#   {
#     ran = subsample_gsgd( xseed = nseed[j], path_sam_table = input_sam_table[i], max_n = 30 )
# 
#     m = match( tab$name[ ran ], fas_id )
#     if( TRUE %in% is.na(m) ){ stop( "mismatch" ) }
# 
#     clade_i  = str_extract( input_sam_table[i], "c[0-9N]{1,6}_" )
#     filename = paste0( "N1TONX/processed/H5_subclade/CN-x2/S4_pH5-", clade_i, "mid", j, ".fasta" )
# 
#     write.fasta( fas_seq[m], fas_id[m], filename )
#   }
#  print(i)
# }

# remove outliers

# toberemoved = read.table( "N1TONX/processed/H5_subclade/CN-x5/outliers.txt",
#                           sep = "\t", stringsAsFactors = FALSE, header = FALSE )
# 
# inputfas    = list.files( "N1TONX/processed/H5_subclade/CN-x4/",
#                           pattern = ".fasta", full.names = TRUE )
# 
# for( f in 1:length(inputfas) )
# {
#   clade = str_match( inputfas[f], "c([0-9N]+)_mid[0-9]" )[,2]
# 
#   rr = which( toberemoved$V1 == clade )
# 
#   if( length(rr) == 0 )
#   {
#     cpname = gsub( ".fasta", "-x5.fasta", inputfas[f] )
#     cmd    = paste0( "cp ", inputfas[f], " ", cpname )
#     system( cmd )
# 
#   }else
#   {
#     seq_id  = fastaEx( inputfas[f] )$id
#     seq_fas = fastaEx( inputfas[f] )$seq
# 
#     m = match( toberemoved$V2[rr], seq_id )
#     m = m[ !is.na(m) ]
# 
#     write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x5.fasta", inputfas[f] ) )
#   }
# }

cat("_BEAST_xml_file_.R")

