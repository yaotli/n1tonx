####### 
#
# This module aims to present the process assigning ecological status to each isolates 
#
# criteria:
#  a. identiy human isolates / major key words (i.e.`wild`)
#  b. assign major domestic birds (i.e. chicken and duck)
#  c. assign other species of animals based on the literature if possible
#  d. examine environment samples based on the literature if possible
#  e. deal with exception 
#
#######

source( "../functions/function.R" )
require( stringr )

# 0 read in ------

fas_CN_x1 = list.files( "N1TONX/processed/H5_subclade/CN-x1", pattern = "^S.*.fasta", full.names = TRUE)

tem_df0 = lapply( as.list(fas_CN_x1),
                 function(x)
                 {
                   names0 = fastaEx( x )$id
                   
                   st1 = gsub( "^[A-Za-z0-9_]+-", "", sapply( strsplit( names0, "__" ), function(x) x[1] ) )
                   st2 = gsub( "A_", "", sapply( strsplit( names0, "__" ), function(x) x[2] ) )
                   
                   clade = str_match( x, "pH5-c([0-9]+)-")[,2]
                   y     = data.frame( ac = st1, name = st2, clade = clade, stringsAsFactors = FALSE )
                   
                   return(y)
                 })

tem_df = do.call( rbind, tem_df0 )
tem_df = tem_df[ order( tem_df$ac ), ]

###### isolate name correlation ######

tem_df$name[ which( tem_df$name == "H5N6" ) ] = c( "Eurasian_Wigeon_Ningxia_476_12_2015", 
                                                   "Ferruginous_Pochard_Ningxia_477_14_2015" )

# 1 cude identification ------

# tem_df$eco  = NA   # D or W
# tem_df$host = NA   # A or M or E
# tem_df$ref = "-"

nonML <- "avian|bird|wildbird|poultry|chicken|duck|dove|pigeon|mallard|goose|environment|water|teal|hawk|eagle|magpie|munia|myna|kestrel|peregrine|crow|sparrow|robin|mesia|gull|egret|swan|shrike|buzzard|heron|quail|pheasant|grebe|starling|swallow|white_eye|cormorant|goldeneye|fowl|pochard|crane|peacock|turkey|falcon|swiftlet|rook|pintail|curlew"

# tem_df$host[ grep( nonML, tem_df$name, ignore.case = TRUE) ] = "A_E"

# not_A_E_i = which( is.na(tem_df$host) )
# not_A_E        = paste0( tem_df$ac, "---", tem_df$name )[ not_A_E_i ] 
# pre_not_A_E    = manual_key( not_A_E )
# pre_not_A_E_df = data.frame( not_A_E, pre_not_A_E, stringsAsFactors = FALSE)

# tem_df$eco[ not_A_E_i ]  = sapply( strsplit( pre_not_A_E_df$pre_not_A_E, ",", fixed = TRUE), function(x) x[1] )
# tem_df$host[ not_A_E_i ] = sapply( strsplit( pre_not_A_E_df$pre_not_A_E, ",", fixed = TRUE), function(x) x[2] )
# tem_df$ref[ not_A_E_i ]  = sapply( strsplit( pre_not_A_E_df$pre_not_A_E, ",", fixed = TRUE), function(x) x[3] )

#
key_wild_i = grep( "wild|migratory", tem_df$name, ignore.case = TRUE )
tem_df$eco[ key_wild_i ]  = "W"
tem_df$host[ key_wild_i ] = "A"

#
key_domestic_i = grep( "domestic", tem_df$name, ignore.case = TRUE )
tem_df$eco[ key_domestic_i ]  = "D"
tem_df$host[ key_domestic_i ] = "A"

# 2 poultry ------

poultry_i = grep( "^duck|^chicken|^goose", tem_df$name, ignore.case = TRUE )
tem_df$eco[ poultry_i ]  = "D"
tem_df$host[ poultry_i ] = "A"

minor_poultry_i = grep( "pigeon|dove|pheasant|peacock|quail|turkey|goose|chicken|duck", tem_df$name, ignore.case = TRUE )
minor_poultry_i = intersect( key_minor_poultry_i,  grep( "W|D", tem_df$eco, invert = TRUE ) )

# minor_poultry        = paste0( tem_df$ac, "---", tem_df$name )[ minor_poultry_i ] 
# pre_minor_poultry    = manual_key( minor_poultry )
# pre_minor_poultry_df = data.frame( minor_poultry, pre_minor_poultry, stringsAsFactors = FALSE )

tem_df$eco[ minor_poultry_i ]  = sapply( strsplit( pre_minor_poultry_df$pre_minor_poultry, ",", fixed = TRUE), function(x) x[1] )
tem_df$host[ minor_poultry_i ] = sapply( strsplit( pre_minor_poultry_df$pre_minor_poultry, ",", fixed = TRUE), function(x) x[2] )
tem_df$ref[ minor_poultry_i ]  = sapply( strsplit( pre_minor_poultry_df$pre_minor_poultry, ",", fixed = TRUE), function(x) x[3] )

# 3 other birds ------

other_birds_i = intersect( grep("W|D|UN", tem_df$eco, invert = TRUE), grep("enviro", tem_df$name, invert = TRUE) )

# other_birds        = paste0( tem_df$ac, "---", tem_df$name )[ other_birds_i ]
# pre_other_birds    = manual_key( other_birds )
# pre_other_birds_df = data.frame( other_birds, pre_other_birds, stringsAsFactors = FALSE )

tem_df$eco[ other_birds_i ]  = sapply( strsplit( pre_other_birds_df$pre_other_birds, ",", fixed = TRUE), function(x) x[1] )
tem_df$host[ other_birds_i ] = sapply( strsplit( pre_other_birds_df$pre_other_birds, ",", fixed = TRUE), function(x) x[2] )
tem_df$ref[ other_birds_i ]  = sapply( strsplit( pre_other_birds_df$pre_other_birds, ",", fixed = TRUE), function(x) x[3] )


# 4 environmemnt ------

environ_i = intersect( grep( "water|enviro", tem_df$name, ignore.case = TRUE ), 
                       grep("W|D|UN", tem_df$eco, invert = TRUE) )

# environ        = paste0( tem_df$ac, "---", tem_df$name )[ environ_i ] 
# pre_environ    = manual_key( environ )
# pre_environ[ 83:86 ] = "D,E,27916476"
# pre_environ_df = data.frame( environ, pre_environ, stringsAsFactors = FALSE  )

tem_df$eco[ environ_i ]  = sapply( strsplit( pre_environ_df$pre_environ, ",", fixed = TRUE), function(x) x[1] )
tem_df$host[ environ_i ] = sapply( strsplit( pre_environ_df$pre_environ, ",", fixed = TRUE), function(x) x[2] )
tem_df$ref[ environ_i ]  = sapply( strsplit( pre_environ_df$pre_environ, ",", fixed = TRUE), function(x) x[3] )

# 5 others ------

eco_i  = grep( "D|W|UN", tem_df$eco, ignore.case = TRUE, invert = TRUE )
host_i = grep( "A|M|E", tem_df$host, ignore.case = TRUE, invert = TRUE )

tem_df$ref[ which( is.na(tem_df$ref) ) ] = "-"
tem_df$host = toupper( tem_df$host )
tem_df$eco  = gsub( "^\\.", "", tem_df$eco )

# UN_i   = grep( "UN", tem_df$eco )
# UN     = paste0( tem_df$ac, "---", tem_df$name )[ UN_i ]
# pre_UN = manual_key( UN )
# UN_ii  = UN_i[ which( pre_UN != "" ) ]

tem_df$eco[ UN_ii ]  = sapply( strsplit( pre_UN[ which( pre_UN != "" ) ], ",", fixed = TRUE), function(x) x[1] )
tem_df$host[ UN_ii ] = sapply( strsplit( pre_UN[ which( pre_UN != "" ) ], ",", fixed = TRUE), function(x) x[2] )
tem_df$ref[ UN_ii ]  = sapply( strsplit( pre_UN[ which( pre_UN != "" ) ], ",", fixed = TRUE), function(x) x[3] )


# write.table( tem_df, "N1TONX/doc/eco/eco_states.tsv",
#              sep = "\t", quote = FALSE, row.names = FALSE )

