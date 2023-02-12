
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

  if iszero(mom1) || iszero(mom2)
    return zero(Basic)
  end # if

  @funs SP

  mom1_array = convert_to_array( mom1 )
  mom2_array = convert_to_array( mom2 )

  result_SP = zero(Basic)
  for pair1 in mom1_array, pair2 in mom2_array
    result_SP += pair1[:num] * pair2[:num] * SP( sort( Basic[pair1[:ki],pair2[:ki]], by=string )... )
  end # for pair1, pair2

  return result_SP

end # function make_SP

###########################################
# Squared momentum 
make_SP( mom::Basic ) = make_SP(mom,mom)
###########################################



###################################################
function make_FV( mom::Basic, rho::Basic )::Basic
###################################################

  @funs FV

  mom_class = SymEngine.get_symengine_class(mom)

  if mom_class == :Symbol

    return FV( mom, rho )

  elseif mom_class == :Mul

    arg_list = get_args(mom)
    @assert length(arg_list) == 2
    first_arg_class = SymEngine.get_symengine_class(arg_list[1])
    the_coeff, the_symbol = first_arg_class == :Integer ? arg_list : reverse(arg_list)
    return the_coeff*FV( the_symbol, rho )

  elseif mom_class == :Add

    arg_list = get_args(mom)
    return sum( map( m_ -> make_FV(m_,rho), arg_list ) )

  else

    error( "Exception found: $(mom), $(rho)" )

  end # if

  return FV(mom,rho)

end # function make_FV


#######################################
function split_SP( expr::Basic )::Basic
#######################################

  @funs FV

  if SymEngine.get_symengine_class( expr ) == :Add
    return (sum ∘ map)( split_SP, get_args(expr) )
  elseif SymEngine.get_symengine_class( expr ) == :Mul
    return (prod ∘ map)( split_SP, get_args(expr) )
  elseif SymEngine.get_symengine_class( expr ) == :Pow
    arg_list = get_args(expr)
    return split_SP( arg_list[1] )^( arg_list[2] )
  elseif SymEngine.get_symengine_class( expr ) == :FunctionSymbol && get_name(expr) == "SP"
    arg_list = get_args(expr)
    return FV(arg_list[1])*FV(arg_list[2])
  else
    return expr
  end # if

  return expr

end # function split_SP




#########################################
function recover_SP( expr::Basic )::Basic
#########################################

  @funs FV, SP

  if SymEngine.get_symengine_class( expr ) == :Add

    return (sum ∘ map)( recover_SP, get_args(expr) )

  elseif SymEngine.get_symengine_class( expr ) == :Mul

    FV_list = filter( x_ -> get_name(x_) == "FV", function_symbols(expr) )
    n_FV = length(FV_list)
    if n_FV == 0
      return expr
    elseif n_FV == 1
      the_FV = first(FV_list)
      mom = (first∘get_args)(the_FV)
      return expr / (the_FV^2) * SP(mom,mom)
    elseif n_FV == 2
      FV1 = FV_list[1]
      mom1 = first(get_args(FV1))
      FV2 = FV_list[2]
      mom2 = first(get_args(FV2))
      return expr / (FV1*FV2) * make_SP(mom1,mom2)
    else
      error( "# of FV is not excepted: "*string(expr) )
    end # if
   
  elseif SymEngine.get_symengine_class( expr ) == :Pow &&
         SymEngine.get_symengine_class( get_args(expr)[1] ) == :FunctionSymbol &&
         get_name( get_args(expr)[1] ) == "FV"

    @assert get_args(expr)[2] == 2
    mom = get_args(get_args(expr)[1])[1]
    return make_SP( mom, mom )

  else

    return expr

  end # if

  return expr

end # function recover_SP




