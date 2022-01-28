# * The score calculation rule (which contains the maximum of three options) is extended with zero as the fourth option.
# * After the dynamic programming table is completed, the maximum score in the entire table is used as the alignment score. The alignment path for this is also identified.
# * The alignment path for the above alignment is "frozen" to zero in the dynamic programming table.
# * A new alignment is calculated, but in the dynamic programming table, the cells that are frozen to zero are not modified. 
#   The resulting alignment is again used to freeze additional cells in the table. This process is repeated until an alignment with a high score (above a predetermined threshold) cannot be found.