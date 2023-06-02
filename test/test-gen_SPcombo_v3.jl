using  Test
using  FeynUtils, SymEngine

@testset "gen_SPcombo_v3" begin
  for q1_rank ∈ 0:3, q2_rank ∈ 0:2, q3_rank ∈ 0:3, num_k ∈ 1:3
    @show [ q1_rank, q2_rank, q3_rank ], [ Basic("k$ii") for ii ∈ 1:num_k ]
    @test gen_SPcombo_v2( [ q1_rank, q2_rank, q3_rank ], [ Basic("k$ii") for ii ∈ 1:num_k ] ) == 
        gen_SPcombo_v3( [ q1_rank, q2_rank, q3_rank ], [ Basic("k$ii") for ii ∈ 1:num_k ] )
  end
end # @testset
