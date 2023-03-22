__precompile__()


module AmpTools

using AbstractAlgebra
using Combinatorics
using Dates
using Groebner
using SHA
using SymEngine


include("BasicMatrix.jl")
include("Type.jl")
include("Conversion.jl")
include("String.jl")
include("Split.jl")
include("Lorentz.jl")
include("Exponent.jl")
include("Groebner.jl")
include("Combo.jl")

export get_Groebner_basis, get_Groebner_basis_v2
export box_message, seq_replace, add_quote, bk_mkdir
export is_FunctionSymbol, is_number, is_class
export to_String_dict, to_Basic_dict, to_Basic, to_String, subs_im
export gen_sorted_str, gen_mma_str
export get_add_vector_noexpand, get_add_vector_expand, get_mul_vector
export get_n_term_noexpand, get_n_term_expand, mul_by_term
export get_n_loop, get_mom_conserv, make_SP, make_FV, split_SP, recover_SP
export get_exponent, get_degree
export split_coeff, drop_coeff, drop_coeff_keep_im
export generate_SPcombo, gen_SPcombo_v2
export iszero_numerical
export get_det, get_adj, get_dot
export calc_sha256


###################
function __init__()
###################
  return nothing
end # function __init__

end # module AmpTools




