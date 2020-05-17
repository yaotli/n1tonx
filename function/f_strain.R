####### 
#
# This module is used to collect eight segments for one particular strain/isolate
#
#######

source( "../functions/function.R" )
require( seqinr )
require( stringr )
require( readxl )
require( purrr )



strain_select = function( tsv_path,
                          fasta_path,
                          genomeset_dat, 
                          na_dat )
{
# # example:
# # folder with sequence metadata of eight segments
# tsv_path   = "/Volumes/EDGE2/db/N1TONX/raw/genomeset/"
# # sequence file for strain selection in GISAID
# fasta_path = "/Volumes/EDGE2/db/N1TONX/raw/genomeset/"
# 
# genomeset_dat = "/Volumes/EDGE2/db/N1TONX/raw/NCBI_ftp/genomeset-20200116.dat"
# na_dat        = "/Volumes/EDGE2/db/N1TONX/raw/NCBI_ftp/influenza_na-20200116.dat"
  
  
  
  input = list.files( tsv_path, 
                      full.names = TRUE,
                      pattern    = "S[0-9].*tsv$" )
  
  seq_files = list.files( fasta_path, 
                          full.names = TRUE,
                          pattern    = ".fasta" )

# 0. read in  #  
  
  df_genome = read.table( genomeset_dat, header = FALSE, sep = "\t", stringsAsFactors = FALSE )
  df_na     = read.delim( na_dat, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE )
  
  df_tsv = do.call( rbind, 
                    lapply( as.list( input ), 
                            function(x) read.table( x, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
                            ) 
                    )
# A. NCBI #
  
  df_tsv_ncbi = df_tsv[ which( df_tsv$isolate_Id == "" ), ]
  
# A.1 match to genomeset.dat #  
  
  # genome_dat may not be accurate
  
  df_tsv_ncbi$isolate_Id = df_genome$V11[ match( df_tsv_ncbi$ac, df_genome$V1 ) ]
  
  ncbi_original_assembled_genome = length(unique(df_tsv_ncbi$isolate_Id))
  ncbi_potentail_adjust_genome   = 0
  ncbi_scattered_genome          = 0
  
# A.2 walk and match #  
  
  unmatched    = which( is.na( df_tsv_ncbi$isolate_Id ) )
  ac_unmatched = df_tsv_ncbi$ac[ unmatched ]
  
  for( i in 1: length(unmatched) )
  {
    if( is.na( df_tsv_ncbi$isolate_Id[ unmatched ][i] ) )
    {
      ac_i = ac_unmatched[i]
      
      s = which( df_na$V1 == ac_i )
      s_u = s + 1000
      s_l = s - 1000
      if( s_u > dim( df_na )[1] ){ s_u = dim( df_na )[1] }
      if( s_l < 1 ){ s_l = 1 }
      
      W = seq(s_l, s_u)
      
      info_pool = paste0( df_na$V4[ W ], df_na$V5[ W ], df_na$V6[ W ],  df_na$V8[ W ] )
      
      # W.i = viruses with similar info
      W.i = W[ which( info_pool == info_pool[ (s-s_l+1) ] ) ]
      
      W.i_ac = df_na$V1[ W.i ]
      
      ac_assigned        = which( !is.na( df_tsv_ncbi$isolate_Id[ match( W.i_ac, df_tsv_ncbi$ac ) ] ) )
      ac_assigned_values = df_tsv_ncbi$isolate_Id[ match( W.i_ac, df_tsv_ncbi$ac ) ][ ac_assigned ]
      
      W.left = W.i[ which( W.i_ac %in% ac_unmatched ) ]
      
      if( (length( ac_assigned ) != 0) & ( length(ac_assigned_values) < length( unique(ac_assigned_values) )*8 ) )
      {
        W.j_ac = df_na$V1[ W.left ]

        df_tsv_ncbi$isolate_Id[ match( W.j_ac, df_tsv_ncbi$ac ) ] =
          unlist( 
            map( W.j_ac, 
                 function(x)
                 {
                   if( is.na( df_tsv_ncbi$isolate_Id[ match( x, df_tsv_ncbi$ac ) ] ) )
                   {
                     y = ac_assigned_values[1]
                     
                   }else{ y = df_tsv_ncbi$isolate_Id[ match( x, df_tsv_ncbi$ac ) ] }
                   
                   return(y)
                 } 
            ) ) 
        ncbi_potentail_adjust_genome = ncbi_potentail_adjust_genome + 1
        
      }else
      {
        if( length( W.left ) > 1 )
        {
          W.j_ac = df_na$V1[ W.left ]
          
          df_tsv_ncbi$isolate_Id[ match( W.j_ac, df_tsv_ncbi$ac ) ] =
            unlist( 
              map( W.j_ac, 
                   function(x)
                   {
                     if( is.na( df_tsv_ncbi$isolate_Id[ match( x, df_tsv_ncbi$ac ) ] ) )
                     {
                       y = paste0( "WM", i ) 
                       
                     }else{ y = df_tsv_ncbi$isolate_Id[ match( x, df_tsv_ncbi$ac ) ] }
                     
                     return(y)
                   } 
              ) )
          ncbi_scattered_genome = ncbi_scattered_genome + 1
          
        }else
        { next() }   
      }
      
print( paste0( i, " among ", length(unmatched) ) )     
      
    }else
    { next() }
    options(warn=2)
    
  }

# A.3 selection #    
  
  multi_segs_isolates0 = 
  sapply( as.list( unique( df_tsv_ncbi$isolate_Id ) ),
          function(x)
          {
            if( is.na(x) )
            {
              return(NA)
            }else
            {
              if( TRUE %in% duplicated( df_tsv_ncbi$seg[ which( df_tsv_ncbi$isolate_Id == x ) ] ) )
              {
                return(x)
              }else
              {
                return(NA)  
              }
            }
          } )
  
  multi_segs_isolates = multi_segs_isolates0[ which( !is.na( multi_segs_isolates0 ) ) ]
  
  for( i in 1: length( multi_segs_isolates ) )
  {
    isolated_i    = which( df_tsv_ncbi$isolate_Id == multi_segs_isolates[i] )
    isolated_ac_i = df_tsv_ncbi$ac[ isolated_i ]
    
    isolated_seg_i = df_na$V3[ match( isolated_ac_i, df_na$V1 ) ]
    isolated_lth_i = df_na$V7[ match( isolated_ac_i, df_na$V1 ) ]
    
    selected_ac = 
    sapply( as.list( seq(1,8) ),
            function(x)
            {
              y = which( isolated_seg_i == x )
                
              if( length( y ) == 1 )
              {
                z = isolated_ac_i[ y ]
                
              }else
              {
                z = isolated_ac_i[ y[ which.max( isolated_lth_i[ y ] ) ] ]
              }
              return(z)
            }
            )
    
    isolated_ac_unselected = match( setdiff( isolated_ac_i, selected_ac ), isolated_ac_i )
    
    df_tsv_ncbi$isolate_Id[ isolated_i[ isolated_ac_unselected ] ] = 
      paste0( "__", df_tsv_ncbi$isolate_Id[ isolated_i[ isolated_ac_unselected ] ], "__" )
  }
  
  print( paste0( "ncbi original assembled genome: ", ncbi_original_assembled_genome ) )
  print( paste0( "ncbi potentail adjust genome: ", ncbi_potentail_adjust_genome ) )
  print( paste0( "ncbi scattered genome: ", ncbi_scattered_genome ) )
  
# B. GISAID #
  
  df_tsv_gisaid = df_tsv[ which( !df_tsv$isolate_Id == "" ), ]
  
  multi_segs_isolates0 = 
    sapply( as.list( unique( df_tsv_gisaid$isolate_Id ) ),
            function(x)
            {
              if( is.na(x) )
              {
                return(NA)
              }else
              {
                if( TRUE %in% duplicated( df_tsv_gisaid$seg[ which( df_tsv_gisaid$isolate_Id == x ) ] ) )
                {
                  return(x)
                }else
                {
                  return(NA)  
                }
              }
            })
  
  multi_segs_isolates2 = multi_segs_isolates0[ which( !is.na( multi_segs_isolates0 ) ) ]
  
  if( length(multi_segs_isolates2) > 0 )
  {
    multi_segs_ls = lapply( as.list( multi_segs_isolates2 ), 
                            function(x)
                            {
                              y = df_tsv_gisaid$ac[ which( df_tsv_gisaid$isolate_Id == x ) ]
                              z = df_tsv_gisaid$seg[ which( df_tsv_gisaid$isolate_Id == x ) ]
                              
                              dup_seg = unique( z[ which( duplicated(z) ) ] )
                              
                              s  = which( z %in% dup_seg )
                              
                              df = data.frame( seg = z[s], ac = y[s], isolate = x, stringsAsFactors = FALSE )
                            } )
    
    multi_segs_ls = do.call( rbind, multi_segs_ls )
    multi_segs_ls$toremove = TRUE
    
    fas_seg_required = sort( unique( multi_segs_ls$seg ) )
    
    for( i in fas_seg_required )
    {
      fas_seq_i = fastaEx( grep( paste0( "S", i, "_" ), seq_files, value = TRUE, ignore.case = TRUE ) )$seq
      fas_id_i  = fastaEx( grep( paste0( "S", i, "_" ), seq_files, value = TRUE, ignore.case = TRUE ) )$id
      
      isolated_i = unique( multi_segs_ls$isolate[ which( multi_segs_ls$seg == i ) ] )
      
      selected_ac = 
      sapply( as.list( isolated_i ),
              function(x)
              {
                 ac_j = which( multi_segs_ls$seg == i & multi_segs_ls$isolate == x ) 
                 j    = match( multi_segs_ls$ac[ac_j], fas_id_i )
                 lth  = sapply(  fas_seq_i[j], 
                                 function(x) length( grep( "[atcg]+", x, ignore.case = TRUE, value = TRUE ) ) )
                 jj   = multi_segs_ls$ac[ac_j][ which.max(lth) ]
              }
              )
      
      multi_segs_ls$toremove[ match( selected_ac, multi_segs_ls$ac  ) ] = FALSE
      
print( paste0( "seg ", i , " done") )
    }
    
    
    m_i = match( multi_segs_ls$ac, df_tsv_gisaid$ac )
    
    df_tsv_gisaid$isolate_Id[m_i][ multi_segs_ls$toremove ] = 
      paste0( "__", df_tsv_gisaid$isolate_Id[m_i][ multi_segs_ls$toremove ], "__" )
  }
  
# C. Export #
  
  out_isolate_id_gisaid = unique(df_tsv_gisaid$isolate_Id)
  out_isolate_id_gisaid = out_isolate_id_gisaid[ which( !startsWith(out_isolate_id_gisaid, "__") ) ]
  
print( paste0("gisaid with ", length(out_isolate_id_gisaid), " isolates") )
  
  gisaid_rows = 
    lapply( as.list( out_isolate_id_gisaid ),
            function(x)
            {
              if( is.na(x) )
              {
                return(x)
              }else
              {
                g_i = which( df_tsv_gisaid$isolate_Id == x )
                
                g_seg = df_tsv_gisaid$seg[g_i]
                g_ac  = df_tsv_gisaid$ac[g_i]
                
                pseudo_df = data.frame( S1 = NA, S2 = NA, S3 = NA, S4 = NA, S5 = NA, S6 = NA, S7 = NA, S8 = NA )
                
                pseudo_df[ ,g_seg ] = g_ac
                
                return(pseudo_df)
              }
            })
  
  out_isolate_id_ncbi = unique(df_tsv_ncbi$isolate_Id)
  out_isolate_id_ncbi = out_isolate_id_ncbi[ which( !startsWith(out_isolate_id_ncbi, "__") ) ]
  
print( paste0("ncbi with ", length(out_isolate_id_ncbi), " isolates") )
  
  ncbi_rows = 
    lapply( as.list( out_isolate_id_ncbi ),
            function(x)
            {
              if( is.na(x) ) 
              {
                return(x)
              }else
              {
                n_i = which( df_tsv_ncbi$isolate_Id == x )
                
                n_seg = df_tsv_ncbi$seg[n_i]
                n_ac  = df_tsv_ncbi$ac[n_i]
                
                pseudo_df = data.frame( S1 = NA, S2 = NA, S3 = NA, S4 = NA, S5 = NA, S6 = NA, S7 = NA, S8 = NA )
                
                pseudo_df[ ,n_seg ] = n_ac
                
                return(pseudo_df)
              }
            })
  
  gisaid_df = do.call( rbind, gisaid_rows )
  ncbi_df   = do.call( rbind, ncbi_rows )
  
  out_df = rbind( gisaid_df, ncbi_df )
  
  rownames( out_df ) = c( out_isolate_id_gisaid, out_isolate_id_ncbi )
  
  filename = paste0( tsv_path, "isolate_ac.tsv" )
  
  write.table( out_df, sep = "\t", file = filename, quote = FALSE )
  
# D. Metadata #

  ls_meta_ncbi = 
  lapply( as.list( seq( 1, dim(ncbi_df)[1] ) ),
          function(x)
          {
            ac_s = match( ncbi_df[x,], df_tsv$ac )  
            ac_s = ac_s[ which( !is.na(ac_s) ) ]
            
            tem_tab  = df_tsv[ ac_s, ][ ,c(2,3,4,5) ]
            tem_comb = paste0( tem_tab$sero, tem_tab$country, tem_tab$year, tem_tab$name ) 
            
            temp     = table( tem_comb )
            con_info = names( temp )[ temp == max(temp) ]
            
            out = tem_tab[ which( tem_comb == con_info )[1], ][, c(1,2,3,4) ]
            return( out )
          } 
          )
  
  df_meta_ncbi = do.call( rbind, ls_meta_ncbi )
  
  
  ls_meta_gisaid = 
  lapply( as.list( out_isolate_id_gisaid ), 
          function(x)
          {
            ac_s = which( df_tsv$isolate_Id == x )[1]
            
            tem_tab = df_tsv[ ac_s, ][ ,c(2,3,4,5) ]
            
            return( tem_tab )
          }
          )
  df_meta_gisaid = do.call( rbind, ls_meta_gisaid )
  
  out_meta_df = rbind( df_meta_gisaid, df_meta_ncbi )
  
  rownames( out_meta_df ) = c( out_isolate_id_gisaid, out_isolate_id_ncbi )
  
  filename2 = paste0( tsv_path, "isolate_meta.tsv" )
  
  write.table( out_meta_df, sep = "\t", file = filename2, quote = FALSE )
  
print("Done")
}



manual_refill_paired_data = function( path_isolate_tsv,
                                      tsv_path )
{

  input = list.files( tsv_path, 
                      full.names = TRUE,
                      pattern    = "S[0-9]_.*.tsv" )
  
  df_tsv = do.call( rbind, 
                    lapply( as.list( input ), 
                            function(x) read.table( x, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
                            ) 
                    )

  ac_tsv = read.table( path_isolate_tsv, header = TRUE, stringsAsFactors = FALSE )
  
  no_non_na = apply( ac_tsv, 1, function(x) length( which( !is.na( x ) ) ) )
  unmatched = names( which( no_non_na != 2) )
  
  ls_seq_ac = 
  lapply( as.list( unmatched ),
          function(x)
          {
            rr   = match( x, rownames(ac_tsv) )
            ac_i = ac_tsv[ rr, ][ ,which( !is.na(ac_tsv[rr,]) ) ]
            
            if( length(ac_i) != 1 | startsWith( x, "WB" ) )
            {
              return( NA )
              
            }else
            {
              y   = which( df_tsv$ac == ac_i )
              y_s = df_tsv$seg[y]
              
              j   = which( df_tsv$name == df_tsv$name[y] )
              j_s = df_tsv$seg[j]
              
              z_s  = setdiff( j_s, y_s )
              z    = which( j_s == z_s )
              z_ac = df_tsv$ac[ j[z] ]
              
              z_ac_e = z_ac[ which(  !z_ac %in% ac_tsv[,z_s] ) ]
              
              if( length( z_ac_e ) == 1 )
              {
                return( c( z_s, z_ac_e ) )
                
              }else
              {
                return( NA )  
              }
            }
          })
        
  for( i in 1: length(unmatched) )
  {
    if( !is.na(ls_seq_ac[[i]][1]) )
    {
      rr = which( rownames(ac_tsv) == unmatched[i] )
      ac_tsv[rr,][, as.numeric(ls_seq_ac[[i]][1]) ] = ls_seq_ac[[i]][2]
      
      print( paste0( "recovered 1 seg:", as.numeric(ls_seq_ac[[i]][1]) ) )
    
    }else
    {
     next()   
    }
  }
  
  if( TRUE %in% duplicated( ac_tsv$S4[ !is.na(ac_tsv$S4) ] ) | TRUE %in% duplicated( ac_tsv$S6[ !is.na(ac_tsv$S6) ] ) ){ stop( "duplicated" ) }
  
  write.table( ac_tsv, gsub( ".tsv", ".2.tsv", path_isolate_tsv ) , sep = "\t", quote = FALSE )
}



update_isolate_ac = function( folder_fasta,
                              path_isolate_tsv )
{
  ac_tsv         = read.table( path_isolate_tsv, header = TRUE, stringsAsFactors = FALSE )
  ac_tsv_updated = ac_tsv
  
  seq_files = list.files( folder_fasta, 
                          full.names = TRUE,
                          pattern    = ".fasta" )
  
  seg_in_folder = as.numeric( str_match( list.files( folder_fasta, pattern = ".fasta" ), "^S([0-9])+")[,2] )
  
  for( i in 1:length(seq_files) )
  {
    fas_id_i = fastaEx( seg_in_folder[i] )$id
    
    ac_tsv_updated[, seg_in_folder[i] ] = 
      ifelse( ac_tsv[, seg_in_folder[i] ] %in% fas_id_i, ac_tsv[, seg_in_folder[i] ], NA ) 
  }
  
  filename = 
  ifelse( grepl( "x[0-9]+.tsv", path_isolate_tsv ), 
          gsub( "x[0-9]+.tsv", paste0( "x", as.numeric( str_match( y, "x([0-9]+).tsv$" )[,2] ), ".tsv" ), path_isolate_tsv),
          gsub( ".tsv", ".x1.tsv", path_isolate_tsv) )
                
          
  write.table( ac_tsv_updated, sep = "\t", file = filename, quote = FALSE )
  
}



#######
# version 20200226
#


