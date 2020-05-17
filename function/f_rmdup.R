####### 
#
# This module is used to remove duplicated sequences considering isolate info at the same time
#
#######

require(seqinr)
require(stringr)
require(lubridate)

fast_rmdup = function( fas_file   = file.choose(),
                       table_file = file.choose() )
{
  
  .test_identical = function( A, B )
  {
    ll = c( nchar(A), nchar(B) )
    
    if( ll[1] > ll[2] )
    {
      b = A 
      a = B 
      
    }else
    {
      a = A
      b = B
    }
    
    l = sort(ll)[1]
    L = sort(ll)[2]
    
    for( k in 1: (L-l+1) )
    { 
      if( substring(b, k, k+l-1) == a )
      {
        return(TRUE)
      }
    }
    return(FALSE)
  }
   
  fas_seq = getSequence( read.fasta( fas_file ) )
  
  fas_id  = attributes( read.fasta( fas_file ) )$names
  fas_id  = gsub( "^[0-9A-Za-z_]+", "", fas_id )
  fas_id  = gsub( "^-", "", fas_id )
  
  tsv_df  = read.table( table_file, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
  
  seq_lth = sapply( fas_seq, 
                    function(x)
                    {
                      s = length( grep( "~|-", x, invert = TRUE, value = TRUE  ) )
                    })
  
  fas_seq_o = fas_seq[ order( seq_lth ) ]
  fas_id_o  = fas_id[ order( seq_lth ) ]
  
  fas_seq_c = lapply( fas_seq_o, 
                      function(x)
                      {
                        y = grep( "~|-", x, invert = TRUE, value = TRUE  ) 
                        y[ !grepl( "[atcg]", y ) ] = "n"
                        
                        return( c2s(y) )
                        #str_replace_all( c2s(y), c("a"="1", "t"="2", "c"="3", "g"="4", "n"="0" ) )
                      } )
  
  mark = rep( 0, length( fas_seq_c ) )
  while( 0 %in% mark )
  {
    g = which( mark == 0 )[1]
    
    sub_g = which( mark == 0 )
    
    print( paste( g, " / ", length(mark) ) )
    
    g.j = which( sapply( fas_seq_c[ sub_g ], function(x){ .test_identical( fas_seq_c[[g]], x ) } ) )
    g.i = sub_g[ g.j ]
    
    if( length(g.i) == 1 & g.i[1] == g )
    {
      mark[ g ] = 3
      
      next()
      
    }else
    {
      
      ii = match( fas_id_o[ g.i ], tsv_df$ac )
      if(  TRUE %in% is.na(ii)  ){ stop( "mismatch" ) }
      
      # length 
      lth_atcg = sapply( fas_seq_o[g.i],
                         function(x)
                         {
                           y = length( grep( "a|t|c|g", x, ignore.case = TRUE, value = TRUE ) )
                         })
      
      if( length( which( lth_atcg == max( lth_atcg ) ) ) < 1 )
      {
        mark[ g.i ] = 1
        mark[ g.i[ which.max( lth_atcg ) ] ] = 3
        
      }else
      {
        df0 = tsv_df[ ii[  lth_atcg == max( lth_atcg )  ], ]
        
        c0 = g.i[ lth_atcg == max( lth_atcg )  ]
        c1 = as.numeric( !grepl( "^EPI", df0$ac ) )*99
        c2 = as.numeric( !df0$country == "UNLOC" )*9
        c3 = nchar( df0$year )
        
        c_sum = c1 + c2 + c3
        
        if( FALSE %in% ( ( max( c_sum ) - c_sum ) == 0 ) )
        {
          df = data.frame( c0, c1, c2, c3 )
          df = df[ order( df$c1, df$c2, df$c3, decreasing = TRUE ), ] 
          
          mark[ g.i ]     = 1
          mark[ df[1,1] ] = 3
          
        }else
        {
          t0 = decimal_date( ymd( df0$year, truncated = 2 ) )
          
          mark[ g.i ]                  = 1
          mark[ g.i[ which.min(t0) ] ] = 3
          
        }
        
      }
      
    }
  }
  
  keep = which( mark == 3 )
  
  write.fasta( sequences = fas_seq_o[ keep ],
               names     = fas_id_o[ keep ],
               file.out  = gsub( ".fasta", "_rmDup.fasta", fas_file ) )
  
  print( paste0( "done; with ", length( which( mark == 1 ) ), " removed" ) )
  
}    


#######
# version 20200222
#


fast_rmdup_v2 = function( fas_file   = file.choose(),
                          table_file = file.choose() )
{
  
  .test_identical = function( A, B )
  {
    ll = c( nchar(A), nchar(B) )
    
    if( ll[1] > ll[2] )
    {
      b = A 
      a = B 
      
    }else
    {
      a = A
      b = B
    }
    
    l = sort(ll)[1]
    L = sort(ll)[2]
    
    for( k in 1: (L-l+1) )
    { 
      if( substring(b, k, k+l-1) == a )
      {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  fas_seq = getSequence( read.fasta( fas_file ) )
  
  fas_id0 = attributes( read.fasta( fas_file ) )$names
  fas_id0 = sapply( str_split( fas_id0, "__"), function(x) x[1] )
  
  fas_id  = gsub( "^[0-9A-Za-z_]+", "", fas_id0 )
  fas_id  = gsub( "^-", "", fas_id )
  
  tsv_df  = read.table( table_file, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
  
  seq_lth = sapply( fas_seq, 
                    function(x)
                    {
                      s = length( grep( "~|-", x, invert = TRUE, value = TRUE  ) )
                    })
  
  fas_seq_o = fas_seq[ order( seq_lth ) ]
  fas_id_o  = fas_id[ order( seq_lth ) ]
  
  fas_seq_c = lapply( fas_seq_o, 
                      function(x)
                      {
                        y = grep( "~|-", x, invert = TRUE, value = TRUE  ) 
                        y[ !grepl( "[atcg]", y ) ] = "n"
                        
                        return( c2s(y) )
                        #str_replace_all( c2s(y), c("a"="1", "t"="2", "c"="3", "g"="4", "n"="0" ) )
                      } )
  
  mark = rep( 0, length( fas_seq_c ) )
  while( 0 %in% mark )
  {
    g = which( mark == 0 )[1]
    
    sub_g = which( mark == 0 )
    
    print( paste( g, " / ", length(mark) ) )
    
    g.j = which( sapply( fas_seq_c[ sub_g ], function(x){ .test_identical( fas_seq_c[[g]], x ) } ) )
    g.i = sub_g[ g.j ]
    
    if( length(g.i) == 1 & g.i[1] == g )
    {
      mark[ g ] = 3
      
      next()
      
    }else
    {
      
      ii = match( fas_id_o[ g.i ], tsv_df$ac )
      if(  TRUE %in% is.na(ii)  ){ stop( "mismatch" ) }
      
      # length 
      lth_atcg = sapply( fas_seq_o[g.i],
                         function(x)
                         {
                           y = length( grep( "a|t|c|g", x, ignore.case = TRUE, value = TRUE ) )
                         })
      
      if( length( which( lth_atcg == max( lth_atcg ) ) ) < 1 )
      {
        mark[ g.i ] = 1
        mark[ g.i[ which.max( lth_atcg ) ] ] = 3
        
      }else
      {
        df0 = tsv_df[ ii[  lth_atcg == max( lth_atcg )  ], ]
        
        c0 = g.i[ lth_atcg == max( lth_atcg )  ]
        c1 = as.numeric( !grepl( "^EPI", df0$ac ) )*99
        c2 = as.numeric( !df0$country == "UNLOC" )*9
        c3 = nchar( df0$year )
        
        c_sum = c1 + c2 + c3
        
        if( FALSE %in% ( ( max( c_sum ) - c_sum ) == 0 ) )
        {
          df = data.frame( c0, c1, c2, c3 )
          df = df[ order( df$c1, df$c2, df$c3, decreasing = TRUE ), ] 
          
          mark[ g.i ]     = 1
          mark[ df[1,1] ] = 3
          
        }else
        {
          t0 = decimal_date( ymd( df0$year, truncated = 2 ) )
          
          mark[ g.i ]                  = 1
          mark[ g.i[ which.min(t0) ] ] = 3
          
        }
        
      }
      
    }
  }
  
  keep = which( mark == 3 )
  
  write.fasta( sequences = fas_seq_o[ keep ],
               names     = fas_id_o[ keep ],
               file.out  = gsub( ".fasta", "_rmDup.fasta", fas_file ) )
  
  print( paste0( "done; with ", length( which( mark == 1 ) ), " removed" ) )
  
}    


fast_rmdup_v3 = function( fas_file   = file.choose(),
                          table_file = file.choose() )
{
  
  .test_identical = function( A, B )
  {
    ll = c( nchar(A), nchar(B) )
    
    if( ll[1] > ll[2] )
    {
      b = A 
      a = B 
      
    }else
    {
      a = A
      b = B
    }
    
    l = sort(ll)[1]
    L = sort(ll)[2]
    
    for( k in 1: (L-l+1) )
    { 
      if( substring(b, k, k+l-1) == a )
      {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  fas_seq = getSequence( read.fasta( fas_file ) )
  
  fas_id0 = attributes( read.fasta( fas_file ) )$names
  
  fas_id_ac = sapply( str_split( fas_id0, "__"), function(x) x[1] )
  
  fas_ac = gsub( "^[0-9A-Za-z_]+", "", fas_id_ac )
  fas_ac = gsub( "^-", "", fas_ac )
  fac_id = gsub( "-[A-Z0-9]+$", "", fas_id_ac )
  
  tsv_df  = read.table( table_file, header = TRUE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE ) 
  
  seq_lth = sapply( fas_seq, 
                    function(x)
                    {
                      s = length( grep( "~|-", x, invert = TRUE, value = TRUE  ) )
                    })
  
  fas_seq_o = fas_seq[ order( seq_lth ) ]
  fas_ac_o  = fas_ac[ order( seq_lth ) ]
  fas_id_o  = fac_id[ order( seq_lth ) ]
  fas_id0_o = fas_id0[ order( seq_lth ) ]
  
  fas_seq_c = lapply( fas_seq_o, 
                      function(x)
                      {
                        y = grep( "~|-", x, invert = TRUE, value = TRUE  ) 
                        y[ !grepl( "[atcg]", y ) ] = "n"
                        
                        return( c2s(y) )
                        #str_replace_all( c2s(y), c("a"="1", "t"="2", "c"="3", "g"="4", "n"="0" ) )
                      } )
  
  mark = rep( 0, length( fas_seq_c ) )
  while( 0 %in% mark )
  {
    g = which( mark == 0 )[1]
    
    sub_g = which( mark == 0 )
    
    print( paste( g, " / ", length(mark) ) )
    
    g.j = which( sapply( fas_seq_c[ sub_g ], function(x){ .test_identical( fas_seq_c[[g]], x ) } ) )
    g.i = sub_g[ g.j ]
    
    if( length(g.i) == 1 & g.i[1] == g )
    {
      mark[ g ] = 3
      
      next()
      
    }else
    {
      
      ii = match( fas_ac_o[ g.i ], tsv_df$ac )
      if(  TRUE %in% is.na(ii)  ){ stop( "mismatch" ) }
      
      # length 
      lth_atcg = sapply( fas_seq_o[g.i],
                         function(x)
                         {
                           y = length( grep( "a|t|c|g", x, ignore.case = TRUE, value = TRUE ) )
                         })
      
      if( length( which( lth_atcg == max( lth_atcg ) ) ) < 1 )
      {
        mark[ g.i ] = 1
        mark[ g.i[ which.max( lth_atcg ) ] ] = 3
        
      }else
      {
        df0    = tsv_df[ ii[  lth_atcg == max( lth_atcg )  ], ]
        df0$id = fas_ac_o[ match( df0$ac, fas_ac_o ) ]
        
        c0 = g.i[ lth_atcg == max( lth_atcg )  ]
        c1 = as.numeric( !grepl( "^EPI", df0$ac ) )*99
        c2 = as.numeric( !grepl( "^UNID", df0$id ) )*99
        #c2 = as.numeric( !df0$country == "UNLOC" )*9
        c3 = nchar( df0$year )
        
        c_sum = c1 + c2 + c3
        
        if( FALSE %in% ( ( max( c_sum ) - c_sum ) == 0 ) )
        {
          df = data.frame( c0, c1, c2, c3 )
          df = df[ order( df$c1, df$c2, df$c3, decreasing = TRUE ), ] 
          
          mark[ g.i ]     = 1
          mark[ df[1,1] ] = 3
          
        }else
        {
          t0 = decimal_date( ymd( df0$year, truncated = 2 ) )
          
          mark[ g.i ]                  = 1
          mark[ g.i[ which.min(t0) ] ] = 3
          
        }
        
      }
      
    }
  }
  
  keep = which( mark == 3 )
  
  write.fasta( sequences = fas_seq_o[ keep ],
               names     = fas_id0_o[ keep ],
               file.out  = gsub( ".fasta", "_rmDup.fasta", fas_file ) )
  
  print( paste0( "done; with ", length( which( mark == 1 ) ), " removed" ) )
  
}    