####### 
#
# This module contains several function to clean the alignments
#
#######

source( "../functions/function.R" )
require( seqinr )
require( stringr )



crude_rm_short_amb = function( folder_raw_fasta )
{
  dir_path = paste0( folder_raw_fasta, "rm_short_amb" )
  
  cmd1 = paste0( "mkdir ", dir_path )
  system( cmd1 )
  cat( paste0( "create a folder at ", dir_path, "\n" ) )
  
  pooledFasta   = list.files( folder_raw_fasta, pattern = "^S[0-9]+.*fasta$", full.names = TRUE )
  pooledFasta_i = list.files( folder_raw_fasta, pattern = "^S[0-9]+.*fasta$" )
  
  if( length( pooledFasta ) > 8 ){ stop( "Include too many fasta files?" ) }
  
  for( p in 1:length(pooledFasta) )
  {
    p.seq = fastaEx( pooledFasta[p] )$seq
    p.id  = fastaEx( pooledFasta[p] )$id
    
    p.lth = sapply( p.seq, length )
    p.med = median( p.lth )/2
    
    p.rm = which( p.lth < p.med )
    
    outseq = p.seq[ - p.rm ]
    outseq = 
      lapply( outseq, 
              function(x)
              {
                y = gsub( "[^acgturywskmdvhb]", "n", x )
                return(y)
              })
    
    filename = gsub( ".fasta", "-rsa.fasta", pooledFasta_i[p] )
    
    write.fasta( sequences = outseq, 
                 names     = p.id[ - p.rm ], 
                 file.out  = paste0( dir_path, "/", filename ) )
    
print( paste0( p, " done with ", length( p.rm ), " seq" ) )
    
  }
print("Done")
}
  


rm_gap_amb <- function( propblank      = 0.95,
                        folder_raw_fasta,
                        exclude        = c(0,0),
                        maxamb         = 1,
                        exclude_except = 4,
                        exclude_ls     = list( HA = c(1000,1200) ) )
{
  dir_path = paste0( folder_raw_fasta, "rm_gap_amb" )
  
  cmd1 = paste0( "mkdir ", dir_path )
  system( cmd1 )
  cat( paste0( "create a folder at ", dir_path, "\n" ) )
  
  fasta   = list.files( folder_raw_fasta, pattern = "^S[0-9]+.*fasta$", full.names = TRUE )
  fasta_i = list.files( folder_raw_fasta, pattern = "^S[0-9]+.*fasta$" )
  
  if( length( fasta ) > 8 ){ stop( "Include too many fasta files?" ) }
  
  for( i in 1:length(fasta) )
  {
    seq_i = fastaEx( fasta[i] )$seq
    id_i  = fastaEx( fasta[i] )$id
    
    amb = which( sapply( seq_i,
                         function(x)
                         {
                           x.l = length( grep( "~|-", x, invert = TRUE, value = TRUE ) )
                           
                           z.s = grep( "a|t|c|g", x, ignore.case = TRUE )[1]
                           z.e = grep( "a|t|c|g", x, ignore.case = TRUE )[ length(grep( "a|t|c|g", x, ignore.case = TRUE )) ]
                           
                           z   = grep( "n", x, ignore.case = TRUE )
                           z.o = length( which( z < z.s ) ) + length( which( z > z.e ) )
                           
                           x.n = length( grep( "~|-|a|t|c|g", x, invert = TRUE, value = TRUE, ignore.case = TRUE ) )
                           return( (x.n-z.o)/x.l > (maxamb/100) )
                         } 
                         ) )
    
    if( length(amb) != 0 )
    {
      seq_j = seq_i[-amb]
      id_j  = id_i[-amb]
      
    }else
    {
      seq_j = seq_i
      id_j  = id_i
    }
    
    seq_matrix = do.call( rbind, seq_j )
    
    fl            = dim( seq_matrix )[1]
    coltoberemove = apply( seq_matrix, 2,
                           function(x)
                           {
                             blank = length( grep( "~|-", x, value = TRUE ))
                            
                             if ( fl*propblank < blank ){ return(1) }else{ return(0) }
                           })
    
    if( i %in% exclude_except ){ exclude = exclude_ls[[ which( exclude_except == i ) ]] }
    
    if( ( exclude[2] > exclude[1] ) & ( exclude[1] > 1 ) ){ coltoberemove[ exclude[1]:exclude[2] ] = 0 }
    
    if( length( which(coltoberemove == 1) ) > 1 ){ cut_matrix = seq_matrix[ ,-which(coltoberemove == 1) ] }else{ cut_matrix = seq_matrix }
    
    seq_cut = as.list( data.frame(t(cut_matrix), stringsAsFactors = FALSE) )
    
    filename = gsub( "_[a-zA-Z]+.fasta$", "-rga.fasta", fasta_i[i] )
 
    write.fasta( sequences = seq_cut, 
                 names     = id_j, 
                 file.out  = paste0( dir_path, "/", filename ) )
print( i )
  }
print("DONE")
}




gapFill <- function( file.dist = file.choose(),
                     e.start   = 1000,
                     e.end     = 1110 )
{
  library(seqinr)

  file      = read.fasta( file.dist )
  seq.name0 = attributes( file )$names
  seq0      = getSequence( file )

  ran      = sample( length( seq.name0 ), 10 )
  seq_test = seq0[ ran ]

  P.codon  =
  sapply( seq_test,
          function(x)
      	  {
           pos = length( grep( "-|~", x[ 1:e.end ], invert = TRUE, value = TRUE ) ) %% 3
      	  })

  df = as.data.frame( table( P.codon ) )
  df = df[ order( df$Freq, decreasing = TRUE ), ]

  if( df$Freq[1] < 0.8 ){ stop("some problem in codon") }

  if( df[1,1] == 1 ){ s.end = e.end - 1 }else if( df[1,1] == 0 ){ s.end = e.end }else{ s.end = e.end + 1 }

  seq <-
    lapply( seq0,
            function(x)
            {
              x1 = grep( "-|~", x[ e.start: s.end ], invert = TRUE, value = TRUE )
              x2 = c( x1, rep( "-", ( s.end - e.start - length(x1) + 1 )) )
              x[ e.start: s.end ] = x2
              return(x)
            }
	    )

  write.fasta( sequences = seq, 
               file.out  = gsub( "\\.fasta|\\.fas", "_gapfill.fasta", file.dist ),
               names     = seq.name0 )
}


#######
# version 2020214
