####### 
#
# This module is used to combine information from two databases 
# and generate a data sheet containing clean 
# assession number, strain name, collection data, serotype and country
#
# Input (in_fas) should include all files to be combine
# 
# naming format
# NCBI
# >{accession}|{strain}|{serotype}|{country}|{year}-{month}-{day}|{segment}
# GISGID
# DNA Accession no.|Isolate name|Type|LOC|Collection date|Segment number
#
#
# [terminal]
# sed -i "" "s/ /_/g" *.fasta
# sed -i "" "s/'//g" *.fasta
#
#######

source( "../functions/function.R" )
require( seqinr )
require( stringr )
require( readxl )
require( purrr )



raw_dat = function( in_fas,
                    in_xls,
                    filecode,
                    name,
                    dir )
{
# #example: 
# in_fas = c( "/Volumes/EDGE2/db/N1TONX/raw/G_p4_H5.fasta",
#             "/Volumes/EDGE2/db/N1TONX/raw/N_p4_H5.fasta" )
# 
# in_xls = c( "/Volumes/EDGE2/db/N1TONX/raw/G_pH5Nx.xls" )
#             
# name = "H5"
# dir    = "./"
# filecode = c( "1G", "1N")
  
  dir_path = paste0( dir, name )
  
  cmd1 = paste0( "mkdir ", dir_path )
  system( cmd1 )
  cat( paste0( "create a folder at ", dir_path, "\n" ) )
  
# A. split according to segment #
  
  for( i in 1: length(in_fas) )
  {
    seq0  = fastaEx( in_fas[i] )$seq
    name0 = fastaEx( in_fas[i] )$id
    
    tailNum        = str_extract( name0, "[0-9]{1,2}$" )
    unique_tailNum = unique( tailNum )
    
    for( j in 1: length( unique_tailNum ) )
    {
      jj = which( tailNum == unique_tailNum[ j ] )
      
      seq_j  = seq0[ jj ]
      name_j = gsub( "\\|[0-9]{1,2}$", "",  name0[ jj ] )
      
      write.fasta( seq_j, name_j, file.out = paste0(  dir_path, "/", "S", unique_tailNum[ j ], "_", filecode[i], 
                                                      ".tem.fasta" ) )
    }
  }
  
  
# B. clean and compile information #
  
  df_xls = do.call( rbind, lapply( as.list( in_xls ), function(x) read_excel( x ) ) )
  df_xls = as.data.frame( df_xls )
  
  n_split = length( unique( str_extract( filecode, "^[0-9]+" ) ) )
  
  tem_files_path = list.files( dir_path, full.names = TRUE )
  tem_files      = list.files( dir_path )
  
  seg   = unique( str_extract( list.files( dir_path ), "S[0-9]+" ) )
  n_seg = length( seg )
  
  gisaid_files = grep( "G.tem.fasta", tem_files )
  
  for( i in 1: n_seg )
  {
    seg_no = as.numeric( str_extract( seg[ i ], "[0-9]+" ) )
    
    S_i        = grep( seg[i], tem_files )
    S_i_gisaid = intersect( S_i, gisaid_files )
    S_i_ncbi   = setdiff( S_i, gisaid_files )
    
# B.1 GISAID #
    
    g_seq0  = do.call( c, map( tem_files_path[ S_i_gisaid ], function(x) fastaEx(x)$seq ) )
    g_name0 = unlist( map( tem_files_path[ S_i_gisaid ], function(x) fastaEx(x)$id ) )
    
    #ac
    g_ac = paste0( "EPI", sapply( str_split( g_name0, "\\|" ), function(x) x[1] ) )
    
    g_tabl_ac = str_extract( df_xls[, seg_no+1 ], "EPI[0-9]+" )
    g_i       = match( g_ac, g_tabl_ac )
    
    na.match = which( is.na( g_i ) )
    if( length( na.match ) > 0 )
    {
      for( m in 1:length( na.match ) )
      {
        g_i[ na.match[m] ] = grep( g_ac[ na.match[m] ], df_xls[, seg_no+1] )
      }
    }
    if( length( which( is.na( g_i ) ) ) > 0 ){ stop("no mathed AC") }
    
    #sero
    g_sero   = str_extract( df_xls$Subtype[ g_i ], "H[0-9]{1,2}N[0-9]{1,2}" )
    na.match = which( is.na(g_sero) )
    if( length( na.match ) > 0 ){ g_sero[ na.match ] = "UNSERO" }
    
    #country
    g_country_0 = sapply( strsplit( df_xls$Location[ g_i ], "/" ), function(x){ x[2] } )
    g_country_0[ grep( ",", g_country_0 ) ] = gsub( ", [A-Za-z ]+", "", g_country_0[ grep( ",", g_country_0 ) ] )
    
    g_country = gsub( " $|^ ", "", g_country_0 )
    g_country = gsub( " ", "_", g_country )
    g_country = gsub( "\\(.*\\)$", "", g_country )
    g_country = gsub( "[[:punct:]]", "_", g_country )
    g_country = gsub( "[_]+", "_", g_country )
    g_country = gsub( "_+$|^_+", "", g_country )
    
    na.match  = which( is.na(g_country) )
    if( length( na.match ) > 0 ){ g_country[ na.match ] = "UNLOC" }
    
    #year
    g_year   = df_xls$Collection_Date[g_i]
    na.match = which( is.na(g_year) )
    if( length( na.match ) > 0 ){ g_year[ na.match ] = "UNYEAR" }
    
    #strain name
    g_name_0 = df_xls$Isolate_Name[g_i]
    
    g_name = gsub( " $|^ ", "", g_name_0 )
    g_name = gsub( " ", "_", g_name )
    g_name = gsub( "\\(.*\\)$", "", g_name )
    g_name = gsub( "[[:punct:]]", "_", g_name )
    g_name = gsub( "[_]+", "_", g_name )
    g_name = gsub( "_+$|^_+", "", g_name )
    
    na.match = which( is.na(g_name) )
    if( length( g_name ) > 0 ){ g_name[ na.match ] = "UNNAME" }
    
    g_df = data.frame( g_ac, g_sero, g_country, g_year, g_name, 
                       g_Id_raw   = g_name0, 
                       isolate_Id = df_xls$Isolate_Id[ g_i ],
                       seg        = seg_no )
    
    
# B.2 NCBI #
    
    n_seq0  = do.call( c, map( tem_files_path[ S_i_ncbi ], function(x) fastaEx(x)$seq ) )
    n_name0 = unlist( map( tem_files_path[ S_i_ncbi ], function(x) fastaEx(x)$id ) )
    
    #ac
    n_Info = strsplit( n_name0, "\\|" ) 
    n_ac   = sapply( n_Info, function(x) x[ 1 ] )
    if( NA %in% grep( "^[A-Z]+[0-9]+$",  n_ac ) ){ stop( "error in NCBI ac." ) }
    
    #sero
    n_sero   = sapply( n_Info, function(x) x[ 3 ] )
    n_sero   = str_extract( n_sero, "[HNhn0-9]{2,6}" )
    
    na.match = which( is.na(n_sero) )
    if( length( na.match ) > 0 ){ n_sero[ na.match ] = "UNSERO" }
    n_sero = toupper( n_sero ) 
    
    #country
    n_country = sapply( n_Info, function(x) x[ 4 ] )
    
    na.match  = which( n_country == "" )
    if( length( na.match ) > 0 ){ n_country[ na.match ] = "UNLOC" }
    
    #year
    n_year = sapply( n_Info, 
                     function(x)
                     {
                       y = x[ 5 ] 
                       z = gsub( "[-]+$", "", y )
                       return(z)
                     })
    
    na.match  = c( which( n_year == "" ), which( is.na( n_year ) ) )
    if( length( na.match ) > 0 ){ n_year[ na.match ] = "UNYEAR" }
    
    #strain name
    n_name_0 = sapply( n_Info, function(x) x[2] )
    
    n_name = gsub( " $|^ ", "", n_name_0 )
    n_name = gsub( " ", "_", n_name )
    n_name = gsub( "\\(.*\\)$", "", n_name )
    n_name = gsub( "[[:punct:]]", "_", n_name )
    n_name = gsub( "[_]+", "_", n_name )
    n_name = gsub( "_+$|^_+", "", n_name )
    
    na.match = which( is.na(n_name) )
    if( length( n_name ) > 0 ){ n_name[ na.match ] = "UNNAME" }
    
    n_df = data.frame( n_ac, n_sero, n_country, n_year, n_name, 
                       n_Id_raw   = n_name0, 
                       isolate_Id = "",
                       seg        = seg_no )
    
    colnames(g_df) = gsub( "g_", "", colnames(g_df) )
    colnames(n_df) = gsub( "n_", "", colnames(n_df) )
    
    df = rbind( g_df, n_df )  
    
    filename = paste0( dir_path, "/S", seg_no, "_", name )
    
    write.table( df, sep = "\t", file = paste0( filename, ".tsv"), quote = FALSE, row.names = FALSE )
    
    write.fasta( sequences = c( g_seq0, n_seq0 ),
                 names     = df$ac,
                 file.out  = paste0( filename, ".fasta" ) )
    
    cat( paste0( "segment ", seg_no, " done", "\n" ) )
  }
  
  system( paste0( "rm ", dir_path, "/*tem*" ) )
}

#######
# version 20200116
#







