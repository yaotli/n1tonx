####### 
#
# This module presents the process of subsampling with H5_subclade sequences
# depending on _gsgd_subsample_.R
# 
#######

# Based on the results of `_gsgd_sumsample_.R`, reference and alternative strinas are 
# determined and colored as previously (purple). Grouped viruses and viruses to be 
# discarded are also colored on the tree. A difference between the followings and the 
# procedures applied to gsgd is here a stand-alone function is used to generate 
# sumsample_tables. 

require( readxl )
require( treeio )
require( ggtree )
require( purrr )
require( stringr )

###### 1 pre-subsampling ######

# sed -i "" "s/,//g" *col2.tre

pre_sumsample = function( path_crude_nex,
                          yr_range     = 3,
                          batchname    = "",
                          col_grouped  = "996633",
                          col_toremove = "ffff00",
                          col_ref      = "ff00ff",
                          col_null     = "000000",
                          ls_outlier   = c("") )
{
  #1 read in 
  
  crude_nex = read.csv( path_crude_nex, stringsAsFactors = FALSE )
  
  taxa.s = grep( x = crude_nex[,1], pattern = "taxlabels" ) + 1
  ntax   = as.numeric( str_match( grep( "ntax", crude_nex[,1], value = TRUE), "(ntax=)([0-9]+)")[,3] )
  taxa.e = taxa.s + ntax - 1
  
  tre_name  = gsub( "\t|\\[.*\\]", "", crude_nex[, 1][taxa.s: taxa.e] )
  tre_col   = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "color=#([0-9a-z]+)" )[,2]
  tre_group = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "name=([0-9]+)" )[,2]
  
  tre_df = data.frame( no = seq( 1, length(tre_name) ), name = tre_name, col = tre_col, group = tre_group, 
                       stringsAsFactors = FALSE )
  
  #2 color 
  
  if( length( ls_outlier ) >1 ){ tre_df$col[ match( outlier$V1, tre_df$name ) ] = col_toremove }
  
  tre_df$group[ tre_df$col == col_toremove ]   = "toRemove"
  tre_df$group[ tre_df$col == col_ref ]        = "ref"
  
  #3 check subgroup
  
  gp_name = unique( tre_df$group )
  gp_name = grep( "^[0-9]", gp_name, value = TRUE )
  
  check_gp_range = map( gp_name, 
                        function(x)
                        {
                          y = which( tre_df$group == x )
                          
                          z = 
                            sapply( str_split( tre_df$name[y], "__" ), 
                                    function(x)
                                    {
                                      i = str_extract(x[5], "^[0-9]{4}")
                                      j = as.numeric( i )
                                      return( j )
                                    })
                          k = ( range(z)[2] - range(z)[1] ) > yr_range
                          
                          return( k )
                        } ) 
  
  check_gp_range = unlist( check_gp_range )               
  
  #
  cat( paste0( "group number : ",  length( gp_name ), "\n" ) )
  cat( paste0( "to be removed : ", length( tre_df$name[ which( tre_df$col == col_toremove ) ] ), "\n" ) )
  cat( paste0( "number of group with range > ", yr_range, " : ", length( which( check_gp_range ) ), "\n" ) )
  cat( paste0( "reference number : ", length( tre_df$name[ which( tre_df$col == col_ref ) ] ), "\n" ) )
  
  #4 sampling year
  tre_df$year = sapply( str_split( tre_df$name, "__" ),
                        function(x)
                        {
                          y = str_extract( x[5], "^[0-9]{4}" )
                          z = as.numeric( y )  
                          return( z )
                        } )
  
  gp_sam_time = unlist( map( gp_name, 
                             function(x)
                             {
                               y   = which( tre_df$group == x )
                               y_i = sort( tre_df$year[y] )
                               y_u = unique( y_i )
                               
                               z = y_u[ which.max( table( y_i ) ) ]
                               return( z )
                             } ) )
  
  tre_df$sam_time = tre_df$year
  for( i in 1: length(gp_sam_time) )
  {
    idx                  = which( tre_df$group == gp_name[i] )
    tre_df$sam_time[idx] = gp_sam_time[i]
  }
  
  tre_df$group[ which( is.na(tre_df$group) ) ] = "n"
  
  unique_year = sort( unique( tre_df$sam_time[ which( tre_df$group != "toRemove" ) ] ) ) 
  
  # 4.5 check lone group
  df_freq_gp = as.data.frame( table( tre_df$group[ tre_df$group %in% gp_name ] ), stringsAsFactors = FALSE )
  lone_gp    = which( df_freq_gp$Freq == 1 )
  if( length(lone_gp) >0 )
  {
    stop( cat( "Group with n = 1 \n ",  tre_df$name[ match(  df_freq_gp$Var1[lone_gp], tre_df$group ) ] ) )
  }
  
  #5 check sample sizes 
  print( unlist(
    map( unique_year,
         function(x)
         {
           x_r   = which( tre_df$sam_time == x & tre_df$group == "ref" )
           n.x_r = length( x_r )
           
           x_n   = which( tre_df$sam_time == x & tre_df$group == "n" )
           n.x_n = length( x_n )
           
           x_g   = tre_df$group[ which( tre_df$sam_time == x & tre_df$group != "n" & tre_df$group != "ref" & tre_df$group != "toRemove" ) ]
           n.x_g = length( unique( x_g ) )
           
           return( paste0( x, ": r = ", n.x_r, "; n = ", n.x_n, "; group = ", n.x_g  ) )
         }) ) )
  
  outname = paste0( gsub( "[-_0-9A-Za-z\\.]+$", "" ,path_crude_nex ), batchname, "_sam_table.tsv" )
  
  write.table( tre_df, sep = "\t", file = outname, quote = FALSE, row.names = FALSE) 
  cat( paste0( "save sample table at: \n", outname ) )
}


###### X1 subsample_gsgd function ###### 

source("../functions/function.R")

subsample_gsgd = function( xseed = 123, 
                           max_n = 25,
                           path_sam_table )
{
  set.seed( xseed )
  
  tre_df = read.table( path_sam_table, sep = "\t", stringsAsFactors = FALSE, header = TRUE )
  
  unique_year = sort( unique( tre_df$sam_time[ which( tre_df$group != "toRemove" ) ] ) ) 
   
  out = 
  unlist(
  map( unique_year,
       function(x)
       {
         x_r   = which( tre_df$sam_time == x & tre_df$group == "ref" )
         n.x_r = length( x_r )
         
         x_n   = which( tre_df$sam_time == x & tre_df$group == "n" )
         n.x_n = length( x_n )
         
         x_g     = which( tre_df$sam_time == x & tre_df$group != "n" & tre_df$group != "ref" & tre_df$group != "toRemove" )
         x_g_u   = unique( tre_df$group[ x_g ] )
         n_x_g_g = length( x_g_u )
         
         if( length(x_g_u) > 0 )
         {
           non_ref = c( x_n, paste0( "g", x_g_u ) )
           
         }else
         {
           non_ref = x_n 
         }
         
         if( ( n.x_r + n.x_n + n_x_g_g ) <  max_n  )
         {
           gp = grep( "^g[0-9]+", non_ref, value = TRUE )
             
           if( length(gp) > 0 )
           {
            sampled_gp = unlist( map( gp, 
                                      function(x)
                                      {
                                        gp_rawname = gsub( "^g", "", x )
                                        
                                        gp_i = which( tre_df$group == gp_rawname )
                                        if( length(gp_i) == 1 )
                                        {
                                          gp_o = gp_i
                                        }else
                                        {
                                          gp_o = sample( gp_i, 1 )   
                                        }
                                        return(gp_o)
                                      }) ) 
            
            sampled = c( x_r, x_n, sampled_gp ) 
            
           }else
           {
             sampled = c( x_r, x_n )    
           }
           
         }else
         {
           n_allowed = max_n - n.x_r
           
           xfactor = sort( seq( 0.1, 0.9, by = 0.1 ), decreasing = TRUE )
            
           n_allowed_0 = length(non_ref)
           j = 1
           while( n_allowed_0 > n_allowed )
           {
             n_allowed_0 = n_allowed_0*xfactor[j]
             j = j+1
           } 
           
           n_allowed_f = ceiling( n_allowed_0 )
           
           smapled0     = sample( non_ref, n_allowed_f )
           sampled_n_gp = c()
           
           if( TRUE %in% grepl( "^g", smapled0 ) )
           {
             sampled0_g = grep( "^g", smapled0, value = TRUE )
             
             sampled_gp = unlist( map( sampled0_g, 
                                       function(x)
                                       {
                                         gp_rawname = gsub( "^g", "", x )
                                         
                                         gp_i = which( tre_df$group == gp_rawname )
                                         
                                         if( length(gp_i) == 1 )
                                         {
                                           gp_o = gp_i
                                         }else
                                         {
                                           gp_o = sample( gp_i, 1 )   
                                         }
                                         return(gp_o)
                                       }) )    
             
             sampled_n_gp = c( sampled_n_gp, sampled_gp )
           }
           
           if( TRUE %in% !grepl( "^g", smapled0 ) )
           {
             sampled0_n   = as.numeric( grep( "^[0-9]+", smapled0, value = TRUE ) ) 
             sampled_n_gp = c( sampled_n_gp, sampled0_n )
             
           }
           
           sampled = c( x_r, sampled_n_gp )
           
         }
         
       }) )
  
  if( TRUE %in% duplicated( out ) ){ stop("dupplicated results") }
  if( TRUE %in% is.na( out ) ){ stop("NA exits") }
  
  return(out)
}
  

###### X2 special sampling scheme ###### 

source("../functions/function.R")

subsample_gsgd_v2 = function( xseed      = 123, 
                              max_n      = 25,
                              path_sam_table,
                              yr_resample = c( 2007,2009 ),
                              yr_prop     = c( 0,0 ) )
{
  set.seed( xseed )
  
  tre_df = read.table( path_sam_table, sep = "\t", stringsAsFactors = FALSE, header = TRUE )
  
  unique_year = sort( unique( tre_df$sam_time[ which( tre_df$group != "toRemove" ) ] ) ) 
  
  out = 
    unlist(
      map( unique_year,
           function(x)
           {
             if( x %in% yr_resample )
             {
               eff_max_n = max_n*yr_prop[ which( yr_resample == x ) ]
               
             }else
             {
               eff_max_n = max_n
             }
             
             if( eff_max_n == 0 ){ return(NULL) }else
             {
               x_r   = which( tre_df$sam_time == x & tre_df$group == "ref" )
               n.x_r = length( x_r )
               
               x_n   = which( tre_df$sam_time == x & tre_df$group == "n" )
               n.x_n = length( x_n )
               
               x_g     = which( tre_df$sam_time == x & tre_df$group != "n" & tre_df$group != "ref" & tre_df$group != "toRemove" )
               x_g_u   = unique( tre_df$group[ x_g ] )
               n_x_g_g = length( x_g_u )
               
               if( length(x_g_u) > 0 )
               {
                 non_ref = c( x_n, paste0( "g", x_g_u ) )
                 
               }else
               {
                 non_ref = x_n 
               }
               
               if( ( n.x_r + n.x_n + n_x_g_g ) <  eff_max_n  )
               {
                 gp = grep( "^g[0-9]+", non_ref, value = TRUE )
                 
                 if( length(gp) > 0 )
                 {
                   sampled_gp = unlist( map( gp, 
                                             function(x)
                                             {
                                               gp_rawname = gsub( "^g", "", x )
                                               
                                               gp_i = which( tre_df$group == gp_rawname )
                                               if( length(gp_i) == 1 )
                                               {
                                                 gp_o = gp_i
                                               }else
                                               {
                                                 gp_o = sample( gp_i, 1 )   
                                               }
                                               return(gp_o)
                                             }) ) 
                   
                   sampled = c( x_r, x_n, sampled_gp ) 
                   
                 }else
                 {
                   sampled = c( x_r, x_n )    
                 }
                 
               }else
               {
                 n_allowed = eff_max_n - n.x_r
                 
                 xfactor = sort( seq( 0.1, 0.9, by = 0.1 ), decreasing = TRUE )
                 
                 n_allowed_0 = length(non_ref)
                 j = 1
                 while( n_allowed_0 > n_allowed )
                 {
                   n_allowed_0 = n_allowed_0*xfactor[j]
                   j = j+1
                 } 
                 
                 n_allowed_f = ceiling( n_allowed_0 )
                 
                 smapled0     = sample( non_ref, n_allowed_f )
                 sampled_n_gp = c()
                 
                 if( TRUE %in% grepl( "^g", smapled0 ) )
                 {
                   sampled0_g = grep( "^g", smapled0, value = TRUE )
                   
                   sampled_gp = unlist( map( sampled0_g, 
                                             function(x)
                                             {
                                               gp_rawname = gsub( "^g", "", x )
                                               
                                               gp_i = which( tre_df$group == gp_rawname )
                                               
                                               if( length(gp_i) == 1 )
                                               {
                                                 gp_o = gp_i
                                               }else
                                               {
                                                 gp_o = sample( gp_i, 1 )   
                                               }
                                               return(gp_o)
                                             }) )    
                   
                   sampled_n_gp = c( sampled_n_gp, sampled_gp )
                 }
                 
                 if( TRUE %in% !grepl( "^g", smapled0 ) )
                 {
                   sampled0_n   = as.numeric( grep( "^[0-9]+", smapled0, value = TRUE ) ) 
                   sampled_n_gp = c( sampled_n_gp, sampled0_n )
                   
                 }
                 
                 sampled = c( x_r, sampled_n_gp )
                 
               }
             }
             
           } ) )
  
  if( TRUE %in% duplicated( out ) ){ stop("dupplicated results") }
  if( TRUE %in% is.na( out ) ){ stop("NA exits") }
  
  return(out)
}


###### 2 subsample ######

# 1 create subsample_tables

# pre_sumsample( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c232-CN_rmDup_col2.tre",
#                batchname = "c232")
# pre_sumsample( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c234-CN_rmDup_col2.tre",
#                batchname = "c234")
# pre_sumsample( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c234N1-CN_rmDup_col2.tre",
#                batchname = "c234N1")
# pre_sumsample( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c2344-CN_rmDup_col2.tre",
#                batchname = "c2344")
# pre_sumsample( "N1TONX/processed/H5_subclade/CN-x1/rmdup/iqtree/S4_pH5-c7-CN_rmDup_col2.tre",
#                batchname = "c7")


# input_sam_table = list.files( "N1TONX/processed/H5_subclade/CN-x2", 
#                               pattern = "sam_table", full.names = TRUE )
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
#     clade_i  = str_extract( input_sam_table[i], "c[0-9N]{1,6}_" )
#     filename = paste0( "N1TONX/processed/H5_subclade/CN-x2/S4_pH5-", clade_i, "sam", j, ".fasta" )
#                        
#     write.fasta( fas_seq[m], fas_id[m], filename ) 
#   }
#  print(i) 
# }









