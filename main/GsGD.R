###### Export for large trees ######

# After sequence level curating, the sequences do not undergo update as genomeset, i.e,
# do not remove sequences that are not paired/indexted in the isolate_table. The 
# alignemnt is then corrected with full-info sequence names. By eyeballing the tree 
# (slc-x1) outliers are identified and removed, resulting in slc-x2 aligment. Isolates 
# in China and Hong Kong are extracted from slc-x2. With function `fast_rmdup_v2`, 
# duplicated sequences in slc-x2-CN are removed. The sequemces name are changed to ac. 
# because `fast_rmdup_v2` is descending from `fast_rmdup` that only takes ac but not 
# isolate ID. Finally, the sequences are subsampled in `_gsgd_subsample_.R`.  

source( "../functions/f_strain.R" )

# 1 copy the files from `raw` to `processed` ------

# cml
# cp -p N1TONX/raw/pH5Nx/align/rm_gap_amb/*_trim.fasta N1TONX/processed/pH5Nx/raw
# for f in $(ls *fasta); do mv "$f" "${f/-*.fasta/-slc.fasta}"; done
# cp -p N1TONX/raw/pH5Nx/isolate_ac.2.tsv N1TONX/processed/pH5Nx/raw
# cp -p N1TONX/raw/pH5Nx/isolate_meta.tsv N1TONX/processed/pH5Nx/raw

# here we do not update the table 
# instead, we simply make a copy - isolate_ac.x1.tsv
# cp isolate_ac.2.tsv isolate_ac.x1.tsv

# 2 update the sequence "name" ------

# slc_fasta   = "N1TONX/processed/pH5Nx//raw/"
# slc_fasta_i = list.files( slc_fasta, full.names = TRUE, pattern = "-slc.fasta" )
# ac_tab_xi   = read.table( paste0( slc_fasta, "isolate_ac.x1.tsv"),
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

# ls_slc_x1 = list.files( "N1TONX/processed/pH5Nx/raw/", pattern = "slc-x1.fasta$", full.names = TRUE )
# ls_tsv    = list.files( "N1TONX/raw/pH5Nx/", pattern = "S[0-9].*tsv$", full.names = TRUE )
# ls_tsv    = c( ls_tsv, rep(ls_tsv[2], 3)  )
# 
# for( r in 1:5)
# {
#   leaf( infile     = ls_slc_x1[r], 
#         tsvfile    = ls_tsv[r], 
#         includeOld = TRUE, 
#         isTree     = FALSE,
#         include    = c( "name", "sero", "country", "year" ))
# }

# for f in $(ls *fasta); do mv "$f" "${f/_rename.fasta/.fasta}"; done
# fasttree

# 3 manual select outliers in each gene and remove ------

# path_slc_x1 = "N1TONX/processed/pH5Nx/slc-x1/"
# path_tsv    = "N1TONX/processed/pH5Nx/slc-x2/outlier.tsv"
# 
# tsv        = read.table( path_tsv, header = FALSE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )
# list_fasta = list.files( path_slc_x1, pattern = ".fasta", full.names = TRUE )
# seg        = as.numeric( str_match( list_fasta, "S([0-9]).*[.fasta]$"  )[,2] )
# 
# for( i in 1: length(list_fasta) )
# {
#   
#   x = which( tsv$V1 == seg[i] )
#   y = tsv$V2[x]
# 
#   fas_seq = fastaEx( list_fasta[i] )$seq
#   fas_id  = fastaEx( list_fasta[i] )$id
# 
#   z = match( y, fas_id )
#   if( length(z) == 0 )
#   {
#     write.fasta( sequences  = fas_seq, names = fas_id,
#                    file.out = gsub( "-x1", "-x2", list_fasta[i] ) )
#   }else
#   {
#     write.fasta( sequences = fas_seq[-z], names = fas_id[-z],
#                  file.out  = gsub( "-x1", "-x2", list_fasta[i] ) )
#   }
# }

# 4 extract CN-HK sequences ------

source( "../functions/function.R" )

# path_tsv    = "N1TONX/raw/pH5Nx/"
# path_slc_x2 = "N1TONX/processed/pH5Nx/slc-x2/"
# 
# input  = list.files( path_tsv, full.names = TRUE, pattern = "^S.*tsv$" )
# df_tsv = do.call( rbind,
#                   lapply( as.list( input ),
#                           function(x) read.table( x, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) ) )
# 
# list_fasta = list.files( path_slc_x2, pattern = ".fasta", full.names = TRUE )
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
#                file.out  = gsub( "-x2", "-CN", list_fasta[i]  ) )
# }


# 5. remove duplicate sequences ------

source( "../functions/f_rmdup.R" )

# fast_rmdup_v2( "N1TONX/processed/pH5Nx/slc-CN/S4_pH5-slc-CN.fasta",
#                "N1TONX/raw/pH5Nx/S4_pH5.tsv" )


# 6. subsample to a big GsGD tree from 1996-2019 with 400-600 HA sequences ------

cat( "../main/_gsgd_subsample_.R" )

# remove two outliers found in TempEst 

# toberemoved = c( "MH341606__A_duck_Shanghai_ZXC03_2017__H5N1__China__2017-11",
#                  "JX878683__A_duck_HuBei_03_2010__H5N5__China__2010-03-28" ) 
# 
# inputfas    = list.files( "N1TONX/processed/pH5Nx/CN-x1/", 
#                           pattern = ".fasta", full.names = TRUE )
# 
# for( f in 1:length(inputfas) )
# {
#   seq_id  = fastaEx( inputfas[f] )$id
#   seq_fas = fastaEx( inputfas[f] )$seq
#   
#   m = match( toberemoved, seq_id)
#   
#   write.fasta( seq_fas[-m], seq_id[-m], gsub( ".fasta", "-x2.fasta", inputfas[f] ) )
# }



