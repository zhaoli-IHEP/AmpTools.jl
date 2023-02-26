__precompile__()


module AmpTools

using Combinatorics
using Dates
using SHA
using SymEngine


include("BasicMatrix.jl")
include("Type.jl")
include("Conversion.jl")
include("String.jl")
include("Split.jl")
include("Lorentz.jl")
include("Exponent.jl")
include("Combo.jl")

export box_message, seq_replace, add_quote, bk_mkdir
export is_FunctionSymbol, is_class
export to_String_dict, to_Basic_dict, to_Basic, to_String
export gen_sorted_str, gen_mma_str
export get_add_vector_noexpand, get_add_vector_expand
export mul_by_term
export convert_to_array
export make_SP, make_FV, split_SP, recover_SP
export get_exponent
export get_degree
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




