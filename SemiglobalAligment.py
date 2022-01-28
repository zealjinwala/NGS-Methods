# Allowing the sequences to have flanking ends without incurring any gap penalties. (aka "free end gaps")

# * First row and first column of the dynamic programming score table is set to zero.
# * After the dynamic programming table is completed, the maximum score in the last row or last column is used as 
#   the alignment score.