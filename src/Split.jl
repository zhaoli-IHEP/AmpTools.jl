#
#
#
#
#
#
####################################################
#"""
#    convert_to_array( 
#        ::Val{:Symbol}, 
#        mom::Basic 
#    )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
#
#This specific conversion is applied only on single symbol momentum, e.g. `k1`.
#
## Examples
#```julia-repl
#julia> using SymEngine, FeAmGen
#
#julia> @vars k1
#(k1,)
#
#julia> FeAmGen.convert_to_array(k1)
#1-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
# (num = 1, ki = k1)
#```
#"""
#function convert_to_array( 
#    ::Val{:Symbol}, 
#    mom::Basic 
#)::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
####################################################
#
#  return NamedTuple{(:num,:ki),Tuple{Basic,Basic}}[ ( num=one(Basic), ki=mom ) ]
#
#end # function convert_to_array
#
#
####################################################
#"""
#    convert_to_array( 
#        ::Val{:Mul}, 
#        mom::Basic 
#    )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
# 
#This specific convertion is applied only on non-unit coefficient with single symbol momentum, e.g. `2*k1`.
# 
## Examples
#```julia-repl
#julia> using SymEngine, FeAmGen
# 
#julia> @vars k1
#(k1,)
# 
#julia> FeAmGen.convert_to_array(2*k1)
#1-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
# (num = 2, ki = k1)
#```
#"""
#function convert_to_array( 
#    ::Val{:Mul}, 
#    mom::Basic 
#)::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
####################################################
# 
#  arg_list = get_args(mom)
#  @assert length(arg_list) == 2
#  first_arg_class = SymEngine.get_symengine_class(arg_list[1])
#  coeff, symbol = first_arg_class == :Integer ? arg_list : reverse(arg_list)
# 
#  return NamedTuple{(:num,:ki),Tuple{Basic,Basic}}[ ( num=coeff, ki=symbol ) ]
# 
#end # function convert_to_array
#
####################################################
#"""
#    convert_to_array( 
#        ::Val{:Add}, 
#         mom::Basic 
#    )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
# 
#This specific convertion is applied only on the summation of non-unit coefficient with single symbol momentum, e.g. `2*k1+3*k2+k3`.
# 
## Examples
#```julia-repl
#julia> using SymEngine, FeAmGen
# 
#julia> @vars k1, k2, k3
#(k1, k2, k3)
# 
#julia> FeAmGen.convert_to_array(2*k1+3*k2+k3)
#3-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
# (num = 1, ki = k3)
# (num = 2, ki = k1)
# (num = 3, ki = k2)
#```
#"""
#function convert_to_array( 
#    ::Val{:Add}, 
#    mom::Basic 
#)::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
####################################################
# 
#  arg_list = get_args(mom)
#  n_solo = length(arg_list)
#  @assert n_solo >= 2
# 
#  return map( arg_list ) do solo_mom
# 
#    solo_arg_list = get_args(solo_mom)
#    @assert length(solo_arg_list) == 2 || is_class(:Symbol,solo_mom) 
# 
#    if length(solo_arg_list) == 2
# 
#      first_arg_class = SymEngine.get_symengine_class(solo_arg_list[1])
#      coeff, symbol = first_arg_class == :Integer ? solo_arg_list : reverse(solo_arg_list)
#      return ( num=coeff, ki=symbol ) # return for map
# 
#    else # solo_mom is :Symbol
# 
#      return ( num=one(Basic), ki=solo_mom ) # return for map
# 
#    end # if
# 
#  end # do solo_mom
# 
# 
#end # function convert_to_array
# 
#
#
####################################################
#"""
#    convert_to_array( mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
#
#This is the interface to different specific convertion functions according to the class of `mom`.
#Now we assume there only there different classes of momentum `mom`. For example, 
#1. `k1`
#2. `2*k1`
#3. `2*k1+3*k2+k3`
#"""
#function convert_to_array( mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
####################################################
#
#  mom_class = SymEngine.get_symengine_class(mom)
#
#  return convert_to_array( Val(mom_class), mom )
#
#end # function convert_to_array

#######################################
@inline function get_n_term_noexpand( 
    expr::Basic 
)::Int64
#######################################

  return is_class(:Add,expr) ? (length∘get_args)(expr) : 1

end # function get_n_term_noexpand

#########################################
@inline function get_n_term_expand( 
    expr::Basic 
)::Int64
#########################################

  expr_expand = expand(expr)
  return is_class(:Add,expr_expand) ? (length∘get_args)(expr_expand) : 1

end # function get_n_term_expand




###########################################
@inline function get_add_vector_noexpand( 
    expr::Basic 
)::Vector{Basic}
###########################################

  return is_class( :Add, expr ) ? get_args(expr) : Basic[ expr ]

end # function get_add_vector_noexpand

#########################################
@inline function get_add_vector_expand( 
    expr::Basic 
)::Vector{Basic} 
#########################################

  return (get_add_vector_noexpand∘expand)(expr)

end # function get_add_vector_expand

###################################
@inline function get_mul_vector( 
    expr::Basic 
)::Vector{Basic}
###################################

  return is_class( :Mul, expr ) ? get_args(expr) : Basic[ expr ]

end # function get_mul_vector



##########################
function mul_by_term(
    expr::Basic,
    mul_fac::Basic
)::Basic
##########################

  term_list = get_add_vector_expand(expr)
  new_result = zero(Basic)
  for one_term in term_list
    new_result += one_term*mul_fac
  end # for one_term

  return new_result

end # function mul_by_term



