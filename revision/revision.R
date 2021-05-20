
source("../functions/function.R")

###### R1 Full-sample dynamics of H5 ######

# 1 copy fas files to the rebuttal/fullsample_dynamics

# 2 extract c234 viruses isolated before 2015

path_c232 = "N1TONX/processed/rebuttal/fullsample_dynamics/S4_pH5-c234-CN_rmDup.fasta"

fas_seq = fastaEx( path_c232 )$seq
fas_id  = fastaEx( path_c232 )$id

info_split    = str_split( string = fas_id, pattern = "__" )
info_split_yr = as.numeric( str_match( lapply( info_split, function(x) x[ length(x) ] ), "[0-9]+" ) )

m = which( info_split_yr < 2015 )
        
# write.fasta( sequences = fas_seq[m], names = fas_id[m],
#              file.out  = gsub( ".fasta$", "-2014.fasta", path_c232 ) )


# 3 extract all CN-rmdup H6-ST2853

path_c6_fas  = "N1TONX/processed/pH6Nx/slc-CN/rmdup/S4_pH6-slc-CN_rmDup.fasta"
path_col_tre = "N1TONX/processed/rebuttal/fullsample_dynamics/H6/S4_pH6-slc-CN_rmDup_colX.tre" 

h6_tre_tab = tagExtra( path_col_tre )

fas_seq = fastaEx( path_c6_fas )$seq
fas_id  = fastaEx( path_c6_fas )$id

id_st2853 = h6_tre_tab$id[ which( h6_tre_tab$tag == "0000ff" ) ]

m = match( id_st2853, fas_id )
which( is.na(m) )

# write.fasta( sequences = fas_seq[m], names = fas_id[m],
#              file.out  = "N1TONX/processed/rebuttal/fullsample_dynamics/H6/S4_pH6-slc-CN_ST2853.fasta" )

# remove one outliers found in TempEst 
# UNID-MH592269__A_goose_China_G2097_2015__H6N2__China__2015-11-07

out = which( fas_id == "UNID-MH592269__A_goose_China_G2097_2015__H6N2__China__2015-11-07" )

m_x = m[ -which(m == out) ]

# write.fasta( sequences = fas_seq[m_x], names = fas_id[m_x],
#              file.out  = "N1TONX/processed/rebuttal/fullsample_dynamics/H6/S4_pH6-slc-CN_ST2853-x1.fasta" )


###### R2 Model comparision ######




