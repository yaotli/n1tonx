####### 
#
# This module facilitates usage of package skygrowth by Erik Volz and Xavier Didelot
# 
#######

# REFERENCE
# Erik Volz and Xavier Didelot, Modeling the growth and decline of pathogen effective 
# population size provides insight into epidemic dynamics and drivers of antimicrobial 
# resistance, Systematic Biology, 2018

# require( devtools )
# install_github("mrc-ide/skygrowth")

growth_df = function( path_anno_trees, 
                      most_recent_time = c( 1900, 2000 ),
                      set_res          = 75,
                      set_tau0         = 0.2, 
                      batch_name       = "c999",
                      names            = c( "A", "B", "C" ),
                      save_file        = FALSE )
{
  require( ggplot2 )
  require( dplyr )
  require( skygrowth )
  require( ape )
  require( ggtree )
  
  t0 = Sys.time()
  
  if( length( most_recent_time ) != length( path_anno_trees ) ){ stop( "no. of most recent time != no. of files"  ) }
  if( length( names ) !=  length( path_anno_trees ) ){ names = paste0( "f", seq( 1, length(path_anno_trees) ) ) }
  
  ls_mcmc = list()
  ls_rate = list()
  
  for( i in 1: length( path_anno_trees ) )
  {
    tr0 = ggtree::read.beast( path_anno_trees[i] )
    tr  = tr0@phylo
    
    #fit    = skygrowth.map( tr, res = set_res, tau0 = set_tau0 )
    mcmcfit = skygrowth.mcmc( tr, res = set_res, tau0 = set_tau0 )
    
    ls_mcmc[[ i ]] = data.frame( mcmcfit$ne_ci, 
                                 Time             = mcmcfit$time + most_recent_time[i], 
                                 Note             = paste0( names[i], "_Ne" ),
                                 Batch            = batch_name, 
                                 stringsAsFactors = FALSE )
    
    ls_rate[[ i ]] = data.frame( mcmcfit$growthrate_ci, 
                                 Time             = mcmcfit$time + most_recent_time[i], 
                                 Note             = paste0( names[i], "_R" ),
                                 Batch            = batch_name, 
                                 stringsAsFactors = FALSE )
  }
  
  df_mcmc = do.call( rbind, ls_mcmc )
  df_rate = do.call( rbind, ls_rate )
  
  colnames( df_mcmc ) = c( "Lower", "Median", "Upper", "Time", "Note", "Batch" )
  colnames( df_rate ) = c( "Lower", "Median", "Upper", "Time", "Note", "Batch" )
  
  f1 = ggplot( ) + geom_line( data = df_mcmc, aes( x = Time, y = Median, color = Note )  ) + theme( legend.position = "none" )
  f2 = ggplot( ) + geom_line( data = df_rate, aes( x = Time, y = Median, color = Note )  ) + theme( legend.position = "none" )
    
  multiplot(f1, f2, ncol = 2)
  
  if( save_file )
  {
    filedir   = gsub( "[-A-Za-z0-9_\\.]+$", "", path_anno_trees[1] )
    filename1 = paste0( filedir, batch_name, "_Ne.tsv" )
    filename2 = paste0( filedir, batch_name, "_rate.tsv" )
    
    write.table( df_mcmc, sep = "\t", file = filename1, quote = FALSE, row.names = FALSE )
    write.table( df_rate, sep = "\t", file = filename2, quote = FALSE, row.names = FALSE )
  }
  
  cat( paste0( "running time: ",  Sys.time() - t0 , " mins") )
  
  outlist          = list( df_mcmc,  df_rate)
  names( outlist ) = c( "Ne", "rate" )
  
  return( outlist )
}




