
###################################################
"""
    make_SP( mom1::Basic, mom2::Basic)::Basic

This function linearly expands the scalar product of `mom1` and `mom2`. 

# Examples
```julia-repl
julia> using SymEngine, FeAmGen

julia> @vars k1, k2, k3
(k1, k2, k3)

julia> FeAmGen.make_SP( 2*k1+3*k2+k3, k1+k2 )
2*SP(k1, k1) + 2*SP(k1, k2) + 3*SP(k2, k1) + 3*SP(k2, k2) + SP(k3, k1) + SP(k3, k2)
```
"""
function make_SP( mom1::Basic, mom2::Basic )::Basic
###################################################

  @funs SP

  mom1_array = convert_to_array( mom1 )
  mom2_array = convert_to_array( mom2 )

  result_SP = zero(Basic)
  for pair1 in mom1_array, pair2 in mom2_array
    result_SP += pair1[:num] * pair2[:num] * SP( sort( Basic[pair1[:ki],pair2[:ki]], by=string )... )
  end # for pair1, pair2

  return result_SP

end # function make_SP

