using SymEngine
using FeynUtils

mat = Basic[ 0 1 2 3; 0 4 5 6 ]
ns_mat = calc_null_space(mat)
@assert all( iszero, mat*ns_mat )

mat = Basic[ 1 1; 2 2; 3 3]
ns_mat = calc_null_space(mat)
@assert all( iszero, mat*ns_mat )

mat = Basic[ 0 0 0; 0 0 0 ]
ns_mat = calc_null_space(mat)
@assert all( iszero, mat*ns_mat )

mat = Basic[1 0; 0 1; 0 0]
ns_mat = calc_null_space(mat)
@assert all( iszero, mat*ns_mat )

mat = Basic[1 1 1; 2 2 2; 3 3 3]
ns_mat = calc_null_space(mat)
@assert all( iszero, mat*ns_mat )

