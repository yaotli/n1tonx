####### 
#
# This module aims to present the reproducible process of subsampling with GsGD sequences
#
# criteria:
#  a. lineage definition strain / vaccine strain
#  b. remove outlier on the temporal regression line
#  c. group sequences with similar geo-host-sero-date info
#  d. random sample the rest in a yearly base
#
#######

require( readxl )
require( treeio )
require( ggtree )
require( purrr )
require( stringr )

###### 1 map the references to the tree ######

# Based on Smith. et al (2015) and the WHO document ( Antigen and genetic characteristics of zoonotic 
# influenza A viruses ), reference/lineage definition strains were identified. After mapping taxa 
# names on the tree (slc-CN-gsgd), unmapped reference strains are manually filled in the ref_strain.tsv
# resulting in ref_strain.x1.tsv for the next step. 


# source of lineage definition strains 
# irv12324-sup-0003-datas2.xlsx
# 201909_zoonotic_vaccinevirusupdate.pdf ------> WHO_vax_ref.tsv

ref_tab = read_excel( "N1TONX/doc/H5_clade_def/irv12324-sup-0003-datas2.xlsx" )

strain0        = gsub( "[[:punct:]]", "_", ref_tab$Strain )
ref_tab$Strain = gsub( "^A_", "", strain0 )

supp_ls = read.table( "N1TONX/doc/H5_clade_def/WHO_vax_ref.tsv", stringsAsFactors = FALSE )
supp_ls = supp_ls[,1]


tre_col0          = read.beast( "N1TONX/processed/pH5Nx/slc-CN/iqtree/S4_pH5-slc-CN-gsgd_col0.tre" )
tre_col0_tiplabel = tre_col0@phylo$tip.label

tre_col0_ac = sapply( str_split( tre_col0_tiplabel, "__" ),
                      function(x)
                      {
                        st0 = x[1]
                        st  = str_extract( st0, "[0-9A-Z]+$" )
                        return(st)
                      } )

tre_col0_name = sapply( str_split( tre_col0_tiplabel, "__" ),
                        function(x)
                        {
                          st0 = x[2]
                          st  = gsub( "A_", "", st0 )
                          return(st)
                        } )

ref_target_name = c( ref_tab$Strain, supp_ls )
ref_target_ac   = c( ref_tab$ID, rep( "XXXXXXX", length(supp_ls) ) )

ref_target_idx = rep( NA, length( ref_target_name ) )
tre_col0_mark  = rep( 0, length( tre_col0_tiplabel ) )

for( i in 1: length( ref_target_name ) )
{
  m1 = match( ref_target_ac[i], tre_col0_ac )

  if( is.na( m1 ) )
  {
    m2 = match( ref_target_name[i], tre_col0_name )

    if( is.na(m2) )
    {
      print( paste( i, ref_target_name[i], sep = " " ) )

    }else if( length(m2) == 1 & !is.na( m2 ) )
    {
      tre_col0_mark[ m2 ]  = 1
      ref_target_idx[i]    = m2

      #print( paste( ref_target_name[i], tre_col0_name[m2], sep = " " ) )
    }

  }else
  {
    tre_col0_mark[m1] = 1
    ref_target_idx[i] = m1

    #print( paste( ref_target_name[i], tre_col0_name[m1], sep = " " ) )
  }
}

tem_df = data.frame( ac = ref_target_ac, name = ref_target_name,
                     m_to_tree = ref_target_idx,  tree_name = tre_col0_tiplabel[ ref_target_idx ] )

# write.table( tem_df, sep = "\t", file = "N1TONX/processed/pH5Nx/slc-CN/iqtree/ref_strain.tsv", 
#              quote = FALSE, row.names = FALSE )

###### 2 color the slc-x2-gdgd (viruses out of China) ######

# To dealt with reference strains out of China, we map taxa names on the slc-x2-gsgd tree. Strain
# "chicken_Viet_Nam_NCVD_15A59_2015" is manually identified, but strains unidentifed are replaced 
# with strains in supp.tsv, which is manually collected to represent close China isolates/clade 
# 2344 viruses/more recent isolates. 


tem_df_x1   = read.table( "N1TONX/processed/pH5Nx/slc-CN/iqtree/ref_strain.x1.tsv", 
                          header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )

slc_x2_gdgd = read.beast( "N1TONX/processed/pH5Nx/slc-x2/iqtree/S4_pH5-slc-x2-gsgd_col0.tre" )

slc_x2_gdgd_tiplabel = slc_x2_gdgd@phylo$tip.label
slc_x2_gdgd_marker   = NA

slc_x2_gdgd_ac = sapply( str_split( slc_x2_gdgd_tiplabel, "__" ),
                         function(x)
                         {
                           st0 = x[1]
                           st  = str_extract( st0, "[0-9A-Z]+$" )
                           return(st)
                         } )
  
slc_x2_gdgd_name = sapply( str_split( slc_x2_gdgd_tiplabel, "__" ),
                           function(x)
                           {
                             st0 = x[2]
                             st  = gsub( "A_", "", st0 )
                             return(st)
                           } )  

for( i in 1: dim( tem_df_x1 )[1] )
{
  if( is.na( tem_df_x1$tree_name[i] ) )
  {
    m1 = match( tem_df_x1$ac[i], slc_x2_gdgd_ac )
    
    if( is.na( m1 ) )
    {
      m2 = match( tem_df_x1$name[i], slc_x2_gdgd_name )
      
      if( is.na(m2) )
      {
        print( paste( i, tem_df_x1$name[i], sep = " " ) )
        
      }else if( length(m2) == 1 & !is.na( m2 ) )
      {
        slc_x2_gdgd_marker[ m2 ]  = 1
      }
      
    }else
    {
      slc_x2_gdgd_marker[ m1 ] = 1
    }
    
  }else{ next() }
}

slc_x2_gdgd@phylo$tip.label[ which( slc_x2_gdgd_marker == 1 ) ] = 
  paste0( slc_x2_gdgd@phylo$tip.label[ which( slc_x2_gdgd_marker == 1 ) ], "[&!color=#ff00ff]" )

# write.beast( slc_x2_gdgd, "N1TONX/processed/pH5Nx/slc-x2/iqtree/S4_pH5-slc-x2-gsgd_col1.tre" )


###### 3 check reference strains in rmdup tree ######

# To narrow down the sampling sequences, rmduped China database (see main/pH5Nx.R) is generated.
# After modifying the taxa name, previously identified references are mapped to the tree. Four 
# alternative strains are selected to replace unmapped (previously defined) references. Lastly,
# CN_rmDup_col0.tre, identical to CN_rmDup_ro_anno2.tre, is colorred with references, resulting 
# in CN_rmDup_col1.tre. 


source("../functions/function.R")

# leaf( infile     = "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_ro.tre",
#       tsvfile    = "N1TONX/raw/pH5Nx/S4_pH5.tsv",
#       includeOld = FALSE, 
#       isTree     = TRUE,
#       include    = c( "ac", "name", "sero", "country", "year"  ) )

# leaf( infile     = "N1TONX/processed/pH5Nx/slc-CN/rmdup/S4_pH5-slc-CN_rmDup.fasta",
#       tsvfile    = "N1TONX/raw/pH5Nx/S4_pH5.tsv",
#       includeOld = FALSE,
#       isTree     = FALSE,
#       include    = c( "ac", "name", "sero", "country", "year"  ) )

# re-append isolate id 

readtre       = read.tree( "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_ro_anno.tre" )
slc_x2_fas_id = fastaEx( "N1TONX/processed/pH5Nx/slc-x2/S4_pH5-slc-x2.fasta")$id

m = match( str_extract( readtre$tip.label, "^[A-Z0-9]+" ),  
           sapply( str_split( slc_x2_fas_id, "__" ), 
                   function(x)
                   {
                     str_extract( x[1], "[A-Z0-9]+$" )
                   }) ) 
   
readtre$tip.label = slc_x2_fas_id[m]

# write.tree( readtre, file = "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_ro_anno2.tre" )                     
      
# check if the references also in rmdup tree

ref1_0 = read.table( "N1TONX/processed/pH5Nx/slc-CN/iqtree/ref_strain.x1.tsv", 
                   header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE )
ref1   = ref1_0$tree_name[ !is.na( ref1_0$tree_name ) ]

ref2_0 = read.table( "N1TONX/processed/pH5Nx/CN-x1/supp.tsv", stringsAsFactors = FALSE )
ref2   = ref2_0$V2

ref = unique( c( ref1, ref2 ) )

rmdup_tree = read.tree( "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_ro_anno2.tre" )

ref[ which( !ref %in% rmdup_tree$tip.label ) ] = 
  c( "WM926-DQ095625__A_Chicken_Yunnan_493_05__H5N1__China__2005",
     "1026572-KM251466__A_duck_Sichuan_NCXJ16_2014__H5N6__China__2014-04-27",
     "1031199-KT963061__A_environment_Yunnan_YN25_2015__H5N6__China__2015-02-16",
     "1026600-KP284957__A_chicken_Shenzhen_433_2013__H5N6__China__2013-12-10" )

# print( which( ! ref %in% rmdup_tree$tip.label ) )
# mark the strains 

cn_rmdup_col0 = read.beast( "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_col0.tre" )

m = match( ref, cn_rmdup_col0@phylo$tip.label)
cn_rmdup_col0@phylo$tip.label[ m ] = 
  paste0( cn_rmdup_col0@phylo$tip.label[ m ], "[&!color=#ff00ff]" )

# write.beast( cn_rmdup_col0, "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_col1.tre" )


###### 4 pre-subsampling ######


# Here CN_rmDup_col1 taxa are manually aggregrated to distinct groups. Each group is colored 
# with brown and labeled with one of the accession number using figtree. Reference strains 
# have been colored with purple in the last step. Yellow color indicates strains to be 
# removed. To read in the tip annotation, periods in the nexus file (CN_rmDup_col2) has to
# be replaced, and the resulting file can no longer be opened in figtree. After reading the 
# tree as dataframe, color codes are translated to text. A few lines of codes are used to 
# examine the data distribution in the tree. To avoid making this section lengthy, we export
# the table as sample_table.tsv, leaving for the subsample function. 

# 1 read in 

# sed -i "" "s/,//g" S4_pH5-slc-CN_rmDup_col3.tre

crude_nex = read.csv( "N1TONX/processed/pH5Nx/slc-CN/rmdup/iqtree/S4_pH5-slc-CN_rmDup_col3.tre",
                        stringsAsFactors = FALSE )

taxa.s = grep( x = crude_nex[,1], pattern = "taxlabels" ) + 1
ntax   = as.numeric( str_match( grep( "ntax", crude_nex[,1], value = TRUE), "(ntax=)([0-9]+)")[,3] )
taxa.e = taxa.s + ntax - 1

tre_name  = gsub( "\t|\\[.*\\]", "", crude_nex[, 1][taxa.s: taxa.e] )
tre_col   = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "color=#([0-9a-z]+)" )[,2]
tre_group = str_match( string  = gsub( "\t", "", crude_nex[, 1][taxa.s: taxa.e] ), pattern = "name=([0-9]+)" )[,2]

tre_df = data.frame( no = seq( 1, length(tre_name) ), name = tre_name, col = tre_col, group = tre_group, 
                     stringsAsFactors = FALSE )


outlier = read.csv( "N1TONX/processed/pH5Nx/CN-x1/outlier.tsv", header = FALSE )

# 2 check subgroup 

gp_name = unique( tre_df$group )[-1] 

check_gp = map( gp_name, 
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
                  k = ( range(z)[2] - range(z)[1] ) > 3
                  
                  return( k )
                } ) 
check_gp = unlist( check_gp )               

#
cat( paste0( "number of group with range > 3 : ", length( which( check_gp ) ), "\n" ) )
cat( paste0( "group number : ",  length( gp_name ), "\n" ) )


# 3 check colors

# 996633 = grouped
# ffff00 = to remove 
# ff00ff = reference 
# 000000 = none 

col_name = unique( tre_df$col )[-1]

tre_df$col[ match( outlier$V1, tre_df$name ) ] = "ffff00"

tre_df$year = sapply( str_split( tre_df$name, "__" ),
                      function(x)
                      {
                        y = str_extract( x[5], "^[0-9]{4}" )
                        z = as.numeric( y )  
                        return( z )
                      } )

#
cat( paste0( "reference number : ",  
             length( tre_df$year[ which( tre_df$col != "ffff00" & tre_df$col == "ff00ff" ) ] ), "\n" ) )

# 4 sample year 

gp_time = unlist( map( gp_name, 
                       function(x)
                       {
                         y   = which( tre_df$group == x )
                         y_i = sort( tre_df$year[y] )
                         y_u = unique( y_i )
                         
                         z = y_u[ which.max( table( y_i ) ) ]
                         return( z )
                       } ) )
               
tre_df$time = tre_df$year
for( i in 1: length(gp_name) )
{
  idx              = which( tre_df$group == gp_name[i] )
  tre_df$time[idx] = gp_time[i]
}

tre_df$group[ tre_df$col == "ffff00" ]       = "toRemove"
tre_df$group[ tre_df$col == "ff00ff" ]       = "ref"
tre_df$group[ which( is.na(tre_df$group) ) ] = "n"

# 5 check sample sizes 

unique_year = sort( unique( tre_df$time[ which( tre_df$group != "toRemove" ) ] ) ) 

# 
unlist(
map( unique_year,
     function(x)
     {
       x_r   = which( tre_df$time == x & tre_df$group == "ref" )
       n.x_r = length( x_r )
       
       x_n   = which( tre_df$time == x & tre_df$group == "n" )
       n.x_n = length( x_n )
       
       x_g   = tre_df$group[ which( tre_df$time == x & tre_df$group != "n" & tre_df$group != "ref" & tre_df$group != "toRemove" ) ]
       n.x_g = length( unique( x_g ) )

       return( paste0( x, ": r = ", n.x_r, "; n = ", n.x_n, "; group = ", n.x_g  ) )
     }) )
     
# write.table( tre_df, sep = "\t", file = "N1TONX/processed/pH5Nx/CN-x1/sample_table.tsv",
#              quote = FALSE, row.names = FALSE)
     
###### X subsample_gsgd function ###### 

source("../functions/function.R")

tre_df = read.table( "N1TONX/processed/pH5Nx/CN-x1/sample_table.tsv", 
                     sep = "\t", stringsAsFactors = FALSE, header = TRUE )

subsample_gsgd = function( xseed = 123, 
                           max_n = 25 )
{
  set.seed( xseed )
  
  unique_year = sort( unique( tre_df$time[ which( tre_df$group != "toRemove" ) ] ) ) 
   
  out = 
  unlist(
  map( unique_year,
       function(x)
       {
         x_r   = which( tre_df$time == x & tre_df$group == "ref" )
         n.x_r = length( x_r )
         
         x_n   = which( tre_df$time == x & tre_df$group == "n" )
         n.x_n = length( x_n )
         
         x_g     = which( tre_df$time == x & tre_df$group != "n" & tre_df$group != "ref" & tre_df$group != "toRemove" )
         x_g_u   = unique( tre_df$group[ x_g ] )
         n_x_g_g = length( x_g_u )
         
         if( length(x_g_u) > 0 )
         {
           non_ref = c( x_n, paste0( "g", x_g_u ) )
         }else
         {
           non_ref = NA   
         }
         
         if( (n.x_r + n.x_n + n_x_g_g) <  max_n  )
         {
           gp = grep( "^g[0-9]+", non_ref, value = TRUE )
             
           if( length(gp) > 0 )
           {
            sampled_gp = unlist( map( gp, 
                                      function(x)
                                      {
                                        gp_rawname = gsub( "^g", "", x )
                                        
                                        gp_i = which( tre_df$group == gp_rawname )
                                        gp_o = sample( gp_i, 1 )
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
                                         gp_o = sample( gp_i, 1 )
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
  
  if( TRUE %in% duplicated( out ) ){ stop() }
  
  return(out)
}
  

###### 5 subsample ###### 

sam1 = subsample_gsgd( xseed = 123 )
sam2 = subsample_gsgd( xseed = 234 )
sam3 = subsample_gsgd( xseed = 456 )
sam4 = subsample_gsgd( xseed = 567 )
sam5 = subsample_gsgd( xseed = 788 )

fas_seq = fastaEx( "N1TONX/processed/pH5Nx/slc-CN/rmdup/S4_pH5-slc-CN_rmDup_rename.fasta" )$seq
fas_id  = fastaEx( "N1TONX/processed/pH5Nx/slc-CN/rmdup/S4_pH5-slc-CN_rmDup_rename.fasta" )$id

for( k in 1:5 )
{
  x = tre_df$name[ sort( get( paste0( "sam", k ) ) ) ]
  x = gsub( "^[A-Z0-9_]+-", "", x )
  m = match( x, fas_id )
  if( TRUE %in%  is.na(m) ){ stop() }else
  {
    outname = gsub( "-slc-CN_rmDup_rename.fasta", paste0( "-gsgd_sam", k, ".fasta" ),
                    "N1TONX/processed/pH5Nx/slc-CN/rmdup/S4_pH5-slc-CN_rmDup_rename.fasta" )
    # write.fasta( fas_seq[m], fas_id[m], outname )  
  }
}










