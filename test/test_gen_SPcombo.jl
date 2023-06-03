using FeynUtils, SymEngine

# Show various results from gen_SPcombo
for q1_rank in 0:3, q2_rank in 0:2, q3_rank in 0:3, num_k in 1:3
  @show [ q1_rank, q2_rank, q3_rank ], [ Basic("k$ii") for ii in 1:num_k ]
  @show gen_SPcombo( [ q1_rank, q2_rank, q3_rank ], [ Basic("k$ii") for ii in 1:num_k ] ) 
end # for

