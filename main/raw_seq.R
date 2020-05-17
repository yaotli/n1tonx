###### Raw samples ######

path_db_raw = "N1TONX/raw/download/"
 
files_db_raw = list.files( path_db_raw )

###### Combine raw data and clean sequence name ######

# The function first splits the sequences according to gene segments. For each gene segment,
# sequence information including accession number, serotype, country, collecion date, and 
# adjusted virus names is extracted from xls tables or sequence names. Sequences and virus 
# metadata are compiled from two databases and the combined data are exported as .fasta file 
# and .tsv file for each segment.

source( "../functions/f_raw_dat.R" )

# 1 genome set

# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(1,3,14,15) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(2,4) ] ),
#          filecode = c( "1G", "2G", "1N", "2N" ),
#          name     = "genomeset",
#          dir      = path_db_raw )
  
# 2 H5 and Nx pair

# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(5,16) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(11) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pH5",
#          dir      = path_db_raw )
# 
# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(10,21) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(11) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pNx",
#          dir      = path_db_raw )

# 3 H6 and Nx pair

# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(6,17) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(12) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pH6",
#          dir      = path_db_raw )
# 
# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(8,19) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(12) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pNx_H6",
#          dir      = path_db_raw )

# 4 H9 and N2 pair

# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(7,18) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(13) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pH9",
#          dir      = path_db_raw )
# 
# raw_dat( in_fas   = paste0( path_db_raw, files_db_raw[ c(9,20) ] ),
#          in_xls   = paste0( path_db_raw, files_db_raw[ c(13) ] ),
#          filecode = c( "1G", "1N" ),
#          name     = "pN2",
#          dir      = path_db_raw )


###### Compile genome data ######

# The function identifies sequences of different segments as one isolate by matching 
# the accession numbers provided by NCBI to genome.dat downloaded from NCBI ftp. The 
# rest of the sequemces with identical information (serotype, country, date, name)
# are also recognized as one isolate and assigned an isolated ID (starts with "WM").
# Sequences already grouped by gemone.dat, if with partial genome (<8), may also be 
# grouped with sequences that had not been assigned by NCBI as the same isolate but 
# with identical information. Sequences of the same segment but assigned with the 
# same isolate ID are selected based on their length. On the other hand, sequences 
# from GISAID are simply grouped by Isolated IDs provided, and then are subjected to
# selection as isolates from NCBI. Metadata of each isolate is determined by consensus 
# of information possessed by eight sequences. Output includes a accession table 
# and a metadata table in .tsv format. 
#
# Another function seaches for isolate IDs with a single sequences in the HA-NA paired 
# data. It then maps the virus name to the whole .tsv file containing sequence name 
# information. (The previous function only searchs identical information nerar by) 
# If the mapped sequences of other segment are ungrouped the isolate_ac table would be 
# updated. 
 
source( "../functions/f_strain.R" )

# # 1 genome set
# 
# strain_select( tsv_path      = "N1TONX/raw/genomeset/",
#                fasta_path    = "N1TONX/raw/genomeset/",
#                genomeset_dat = "N1TONX/raw/NCBI_ftp/genomeset-20200116.dat",
#                na_dat        = "N1TONX/raw/NCBI_ftp/influenza_na-20200116.dat" )
# 
# # 2 H5 and Nx pair
# 
# # cml
# # mkdir pH5Nx pH6N1 pH9N2 
# # cp pH5/* pH5Nx
# # cp pNx/* pH5Nx 
# # ...
# 
# strain_select( tsv_path      = "N1TONX/raw/pH5Nx/",
#                fasta_path    = "N1TONX/raw/pH5Nx/",
#                genomeset_dat = "N1TONX/raw/NCBI_ftp/genomeset-20200116.dat",
#                na_dat        = "N1TONX/raw/NCBI_ftp/influenza_na-20200116.dat" )
# 
# manual_refill_paired_data( path_isolate_tsv = "N1TONX/raw/pH5Nx/isolate_ac.tsv",
#                            tsv_path         = "N1TONX/raw/pH5Nx/" )
# 
# 
# 3 H6 and N1 pair
# 
# strain_select( tsv_path      = "N1TONX/raw/pH6Nx/",
#                fasta_path    = "N1TONX/raw/pH6Nx/",
#                genomeset_dat = "N1TONX/raw/NCBI_ftp/genomeset-20200116.dat",
#                na_dat        = "N1TONX/raw/NCBI_ftp/influenza_na-20200116.dat" )
# 
# manual_refill_paired_data( path_isolate_tsv = "N1TONX/raw/pH6Nx/isolate_ac.tsv",
#                            tsv_path         = "N1TONX/raw/pH6Nx/" )
# 
# # 4 H9 and N2 pair
# 
# strain_select( tsv_path      = "N1TONX/raw/pH9N2/",
#                fasta_path    = "N1TONX/raw/pH9N2/",
#                genomeset_dat = "N1TONX/raw/NCBI_ftp/genomeset-20200116.dat",
#                na_dat        = "N1TONX/raw/NCBI_ftp/influenza_na-20200116.dat" )
# 
# manual_refill_paired_data( path_isolate_tsv = "N1TONX/raw/pH9N2//isolate_ac.tsv",
#                            tsv_path         = "N1TONX/raw/pH9N2/" )


###### Sequence level cleaning for genomeset ######

# The `crude_rm_short_amb` function removes sequences with length shorter than median value 
# of the sequence lengths in the fasta file. The function also replace non-IUPAC nucleotide 
# symbols with 'n'. Then H5 and N1 sequences are isolated from HA and NA files before 
# alignment. The `rm_gap_amb` eliminates gaps generated by outliers and remove sequences 
# with ambiguous nucleotide exceeding 5%. HA alignment is subjected to function `gapFill`
# to maintain the sequences near cleavage site. Alignments are then manually trimmed to 
# coding region and sequences shorther than 0.9 alignment length are filtered out. The 
# subsequent sequence ID are used to extracted raw sequences to prepare a new file. New
# set of sequences are finally aligned again and timmed with `rm_gap_amb` and `gapFill`
# as previously.

source( "../functions/f_seq_pruing.R" )

# crude_rm_short_amb( folder_raw_fasta = "N1TONX/raw/genomeset/" ) 
#
# # split H5 and N1 out from the pool 
# 
# path_fas = c( "N1TONX/raw/genomeset/rm_short_amb/S4_genomeset-rsa.fasta",
#               "N1TONX/raw/genomeset/rm_short_amb/S6_genomeset-rsa.fasta" )
# path_tsv = c( "N1TONX/raw/genomeset/S4_genomeset.tsv",
#               "N1TONX/raw/genomeset/S6_genomeset.tsv" )
# outname  = c( "N1TONX/raw/genomeset/rm_short_amb/S4_genomeset_H5-rsa.fasta",
#               "N1TONX/raw/genomeset/rm_short_amb/S6_genomeset_N1-rsa.fasta" )
# 
# sero_x   = c( "H5", "N1" )
# 
# for( i in 1:2 )
# {
#   tsv     = read.table( path_tsv[i], header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
#   
#   fas_seq = fastaEx( path_fas[i] )$seq
#   fas_id  = fastaEx( path_fas[i] )$id
#   
#   m   = match( fas_id, tsv$ac )
#   m_i = grep( sero_x[i], tsv$sero[m] )
#   
#   write.fasta( fas_seq[m_i], fas_id[m_i], outname[i]  )
# }

## ./shell/mafft
# rm_gap_amb( propblank        = 0.95,
#             folder_raw_fasta = "N1TONX/raw/genomeset/rm_short_amb_align/",
#             exclude          = c(0,0),
#             maxamb           = 5,
#             exclude_except   = 4,
#             exclude_ls       = list( HA = c(1079,1145) ) )
# gapFill( file.dist = "N1TONX/raw/genomeset/rm_short_amb_align/rm_gap_amb/S4_genomeset_H5-rsa-rga.fasta",
#          e.start   = 1041, 
#          e.end     = 1107 )
# 
# ### 
# #
# # Manual trim and inspect sequence
# #
# 
# folder_trimmed_fasta = "N1TONX/raw/genomeset/rm_short_amb_align/rm_gap_amb/"
# folder_fasta         = "N1TONX/raw/genomeset/rm_short_amb/"
# 
# trimmed_fasta = list.files( folder_trimmed_fasta, pattern = "_trim.fasta", full.names = TRUE )
# fasta         = list.files( folder_fasta, pattern = ".fasta", full.names = TRUE )
# 
# # sequence length filtering
# for( i in 1:8 )
# {
#   fasta_t   = trimmed_fasta[i]
#   fasta_raw = grep( paste0( "S", i, "_" ), fasta, value = TRUE )[1]
#   
#   t_id  = fastaEx( fasta_t )$id
#   t_seq = fastaEx( fasta_t )$seq
#   
#   t_w   = unique( sapply( t_seq, length ) )
#   t_lth = sapply( t_seq,
#                   function(x)
#                   {
#                     y = grep( "~|-", x, invert = TRUE, value = TRUE )
#                     z = length(y)
#                     return( z )
#                   })
#                     
#   tokeep = which( t_lth >= t_w*0.90 )
# 
#   raw_id  = fastaEx( fasta_raw )$id
#   raw_seq = fastaEx( fasta_raw )$seq
#   
#   m_i = match( t_id[tokeep], raw_id )
#   
#   if( length( which( is.na( m_i ) ) ) > 0 ){ stop( "mismatch" ) } 
#   
#   out_id  = t_id[tokeep]
#   out_seq = raw_seq[ m_i ]
#   
#   filename = gsub( "_trim.fasta$", "_tobeReAlign.fasta", fasta_t )
#   
#   write.fasta( sequences = out_seq,
#                names     = out_id, 
#                file.out  = filename )
#   
#   segment_lth_note = paste0( "echo ", t_w, " >> ", "seg_lth.txt"  )
#   system( segment_lth_note )
# }
# 
# 
# # ./shell/mafft
# 
# # trim again 
# 
# rm_gap_amb( propblank        = 0.95,
#             folder_raw_fasta = "N1TONX/raw/genomeset/align/",
#             exclude          = c(0,0),
#             maxamb           = 50,                                              #<------------ note here
#             exclude_except   = 4,
#             exclude_ls       = list( HA = c(1078,1147) ) )
# gapFill( file.dist = "N1TONX/raw/genomeset/align/rm_gap_amb/S4_genomeset_H5-rsa-rga-rga.fasta",
#          e.start   = 1041, 
#          e.end     = 1110 )
# 
# ### 
# #
# # Manual trim and inspect sequence
# #

###### Sequence level cleaning for paired data ######

# Sequence quality filtering is conducted as `genomeset` dataset described in step 3.
# For pH5Nx dataset, NA alignment is split to N1, N2, N6, N8, and for pH6Nx datasets
# NA alignment is split to N1, N2, N5, N6, N8. 

source( "../functions/f_seq_pruing.R" )

# 1 remove short sequences 

# crude_rm_short_amb( folder_raw_fasta = "N1TONX/raw/pH5Nx/" )
# crude_rm_short_amb( folder_raw_fasta = "N1TONX/raw/pH6Nx/" ) 
# crude_rm_short_amb( folder_raw_fasta = "N1TONX/raw/pH9N2/" ) 

# 2 split subtypes ( H5 and H6 )

# path_fas = c( "N1TONX/raw/pH5Nx/rm_short_amb/S6_pNx-rsa.fasta",
#               "N1TONX/raw/pH6Nx/rm_short_amb/S6_pNx_H6-rsa.fasta" )
# path_tsv = c( "N1TONX/raw/pH5Nx/S6_pNx.tsv",
#               "N1TONX/raw/pH6Nx/S6_pNx_H6.tsv" )
# 
# for( i in 1:2 )
# {
#   tsv     = read.table( path_tsv[i], header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )
# 
#   fas_seq = fastaEx( path_fas[i] )$seq
#   fas_id  = fastaEx( path_fas[i] )$id
# 
#   n_NA = unique( tsv$sero )
#   i_NA = gsub( "^H5|H6", "", n_NA ) 
#   
#   for( j in 1: length(n_NA) )
#   {
#     
#     m0 = match( tsv$ac[ which( tsv$sero == n_NA[j] ) ], fas_id )
#     m  = m0[ !is.na( m0 ) ]
#     
#     if( length(m) < 10 ){ next() }else
#     {
#       write.fasta( fas_seq[m], fas_id[m], gsub(  "pNx", paste0( "p", i_NA[j]), path_fas[i] )  ) 
#     } 
#   }
# }

# 3 align and trim 

# ./shell/mafft

# #H5
# rm_gap_amb( propblank        = 0.95,
#             folder_raw_fasta = "N1TONX/raw/pH5Nx/rm_short_amb_align/",
#             exclude          = c(0,0),
#             maxamb           = 5,
#             exclude_except   = 1,
#             exclude_ls       = list( HA = c(1155,1238) ) )
# 
# gapFill( file.dist = "N1TONX/raw/pH5Nx/rm_short_amb_align/rm_gap_amb/S4_pH5-rsa-rga.fasta",
#          e.start   = 1041,
#          e.end     = 1122 )
# 
# #H6
# rm_gap_amb( propblank        = 0.95,
#             folder_raw_fasta = "N1TONX/raw/pH6Nx/rm_short_amb_align/",
#             exclude          = c(0,0),
#             maxamb           = 5,
#             exclude_except   = "" )
# 
# #H9
# rm_gap_amb( propblank        = 0.95,
#             folder_raw_fasta = "N1TONX/raw/pH9N2/rm_short_amb_align/",
#             exclude          = c(0,0),
#             maxamb           = 5,
#             exclude_except   = "" )
###
#
# Manual trim and inspect sequence
#

# 4 filter out short sequence with cutting value - 90%

# folder_trimmed_fasta = c( "N1TONX/raw/pH5Nx/rm_short_amb_align/rm_gap_amb/",
#                           "N1TONX/raw/pH6Nx/rm_short_amb_align/rm_gap_amb/",
#                           "N1TONX/raw/pH9N2/rm_short_amb_align/rm_gap_amb/" )
# 
# folder_fasta         = c( "N1TONX/raw/pH5Nx/",
#                           "N1TONX/raw/pH6Nx/",
#                           "N1TONX/raw/pH9N2/" )
# 
# for ( i in 1:3 )
# {
#   trimmed_fasta = list.files( folder_trimmed_fasta[i], pattern = "_trim.fasta", full.names = TRUE )
#   fasta         = list.files( folder_fasta[i], pattern = ".fasta", full.names = TRUE )
#   
#   for( j in 1: length(trimmed_fasta) )
#   {
#     fasta_t   = trimmed_fasta[j]
#     
#     nseg      = str_extract( trimmed_fasta[j], "S[0-9]" )
#     fasta_raw = grep( nseg, fasta, value = TRUE )[1]
#     
#     t_id  = fastaEx( fasta_t )$id
#     t_seq = fastaEx( fasta_t )$seq
#     
#     t_w   = unique( sapply( t_seq, length ) )
#     t_lth = sapply( t_seq,
#                     function(x)
#                     {
#                       y = grep( "~|-", x, invert = TRUE, value = TRUE )
#                       z = length(y)
#                       return( z )
#                     })
#     
#     tokeep = which( t_lth >= t_w*0.90 )
#     
#     raw_id  = fastaEx( fasta_raw )$id
#     raw_seq = fastaEx( fasta_raw )$seq
#     
#     m_i = match( t_id[tokeep], raw_id )
#     
#     if( length( which( is.na( m_i ) ) ) > 0 ){ stop( "mismatch" ) }
#     
#     out_id  = t_id[tokeep]
#     out_seq = raw_seq[ m_i ]
#     
#     filename = gsub( "_trim.fasta$", "_tobeReAlign.fasta", fasta_t )
#     
#     write.fasta( sequences = out_seq,
#                  names     = out_id,
#                  file.out  = filename ) 
#   }
#   print( i )
# }

# 3 re-align and trim 

# ./shell/mafft

#H5
rm_gap_amb( propblank        = 0.99,
            folder_raw_fasta = "N1TONX/raw/pH5Nx/align/",
            exclude          = c(0,0),
            maxamb           = 50,
            exclude_except   = 1,
            exclude_ls       = list( HA = c(1108, 1204) ) )

gapFill( file.dist = "N1TONX/raw/pH5Nx/align/rm_gap_amb/S4_pH5-rsa-rga-rga.fasta",
         e.start   = 1041,
         e.end     = 1135 )

#6
rm_gap_amb( propblank        = 0.99,
            folder_raw_fasta = "N1TONX/raw/pH6Nx/align/",
            exclude          = c(0,0),
            maxamb           = 50,
            exclude_except   = "" )

#H9
rm_gap_amb( propblank        = 0.99,
            folder_raw_fasta = "N1TONX/raw/pH9N2/align/",
            exclude          = c(0,0),
            maxamb           = 50,
            exclude_except   = "" )

###
#
# Manual trim and inspect sequence
#

