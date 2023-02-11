__precompile__()


module SymEngineExt

using SymEngine


include("BasicMatrix.jl")
include("Type.jl")
include("Conversion.jl")
include("String.jl")
include("Split.jl")
include("Lorentz.jl")
include("Exponent.jl")
include("Combo.jl")

export box_message, sequential_replace
export is_FunctionSymbol
export to_String_dict, to_Basic_dict
export gen_sorted_str, gen_mma_str
export get_add_vector_noexpand, get_add_vector_expand
export mul_by_term
export convert_to_array
export make_SP, make_FV, split_SP, recover_SP
export get_exponent
export get_degree
export drop_coeff, drop_coeff_keep_im
export generate_SPcombo, gen_SPcombo_v2
export iszero_numerical
export get_det, get_adj, get_dot






###################
function __init__()
###################
  return nothing
end # function __init__

end # module SymEngineExt




