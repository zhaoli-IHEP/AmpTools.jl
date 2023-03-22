

###########################
function get_n_loop(
  loop_den_list::Vector{Basic},
)::Int64
###########################

  @assert all( x->is_FunctionSymbol(x)&&get_name(x)=="Den", loop_den_list )

  mom_list = map( first∘get_args, loop_den_list )
  n_loop = (length∘unique∘filter)( x->(first∘string)(x)=='q', free_symbols(mom_list) )
  return n_loop

end # function get_n_loop


##################################
function get_mom_conserv(
    n_inc::Int64,
    ext_mom_list::Vector{Basic}
)::Pair{Basic,Basic}
##################################

  if iszero(n_inc)
    # In the scalar integral, mom_conserv is supposed to have been implemented.
    @vars k1
    return k1 => k1
  end # if

  # momentum conservation
  Kn_expr = expand( 2*sum(ext_mom_list[1:n_inc])-sum(ext_mom_list)+ext_mom_list[end] )
  mom_conserv = ext_mom_list[end] => Kn_expr
  return mom_conserv

end # function get_mom_conserv




######################
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
function make_SP(
    mom1::Basic, 
    mom2::Basic 
)::Basic 
######################

  if iszero(mom1) || iszero(mom2)
    return zero(Basic)
  end # if

  @funs SP
  term_list = get_add_vector_expand(mom1*mom2)
  result_sp = zero(Basic)
  for one_term in term_list
    var_list = sort( free_symbols(one_term), by=string )
    @assert length(var_list) in [1,2]
    dict = Dict( var_list .=> one(Basic) )
    the_coeff = subs( one_term, dict )
    the_sp = length(var_list) == 1 ? SP(var_list...,var_list...) : SP( var_list... )
    result_sp += the_coeff*the_sp
  end # for one_term
  
  return result_sp

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
    return sum( map( x -> make_FV(x,rho), arg_list ) )

  else

    error( "Exception found: $(mom), $(rho)" )

  end # if

  return FV(mom,rho)

end # function make_FV


#######################################
# Split SP(k1,k2) => FV(k1)*FV(k2)
function split_SP( expr::Basic )::Basic
#######################################

  @funs FV

  if is_class( :Add, expr )

    return (sum ∘ map)( split_SP, get_args(expr) )

  elseif is_class( :Mul, expr ) 

    return (prod ∘ map)( split_SP, get_args(expr) )

  elseif is_class( :Pow, expr ) 

    base, xpt = get_args(expr)
    return split_SP(base)^xpt

  elseif is_FunctionSymbol(expr) && get_name(expr) == "SP"

    mom1, mom2 = get_args(expr)
    return FV(mom1)*FV(mom2)

  else

    return expr

  end # if

  return expr

end # function split_SP




#########################################
function recover_SP( expr::Basic )::Basic
#########################################

  @funs FV, SP

  if is_class( :Add, expr )

    return (sum ∘ map)( recover_SP, get_args(expr) )

  elseif is_class( :Mul, expr )

    FV_list = filter( x -> get_name(x) == "FV", function_symbols(expr) )
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
      error( "# of FV is not excepted: $(expr)" )
    end # if
   
  elseif is_class( :Pow, expr ) &&
         is_FunctionSymbol( get_args(expr)[1] ) &&
         (get_name∘first∘get_args)(expr) == "FV"

    @assert (last∘get_args)(expr) == 2
    mom = (first∘get_args∘first∘get_args)(expr)
    return make_SP( mom, mom )

  else

    return expr

  end # if

  return expr

end # function recover_SP


