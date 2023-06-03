

#############################
@inline function is_FunctionSymbol(
    expr::Basic
)::Bool
#############################
  return SymEngine.get_symengine_class(expr) == :FunctionSymbol 
end # function is_FunctionSymbol

#############################
@inline function is_number(
    expr::Basic
)::Bool
#############################
  return SymEngine.get_symengine_class(expr) in [:Integer,:Rational,:Complex]
end # function is_number



#############################
@inline function is_class(
    class::Symbol,
    expr::Basic
)::Bool
#############################
  return SymEngine.get_symengine_class(expr) == class
end # function is_class




############################
function iszero_numerical(
    expr::Basic
)::Bool
############################

  if SymEngine.get_symengine_class(expr) in [:Integer, :Rational, :Complex]
    return iszero(expr)
  end # if

  var_list = free_symbols( expr )
  @assert !isempty(var_list)
  n_var = length(var_list)

  eval_dict = Dict( var_list .=> zeros(Int64,n_var) )

  n_repeat = 4^length(var_list)
  for repeat in 1:n_repeat
    eval_dict = Dict( var_list .=> map( x->mod(x,256)+16, rand(Int64,n_var) ) )

    try
      eval_expr = subs( expr, eval_dict )

      if isnan(eval_expr) || findfirst( "zoo", string(eval_expr) ) != nothing
        printstyled( ">> Found exception for numerical evalulation >>\n", color=:light_yellow )
        eval_dict = Dict( var_list .=> map( x->mod(x,256)+16, rand(Int64,n_var) ) )
        eval_expr = subs( expr, eval_dict )
        @assert !isnan(eval_expr) && findfirst( "zoo", string(eval_expr) ) == nothing
      end # if

      if !iszero(eval_expr)
        return false
      end # if

    catch err
      printstyled( ">> Found error from subs in iszero_numerical, now fix it. >>\n", color=:light_yellow )
      eval_dict = Dict( var_list .=> map( x->mod(x,256)+16, rand(Int64,n_var) ) )
      eval_expr = subs( expr, eval_dict )
      if !iszero(eval_expr)
        return false
      end # if

    end # try

  end # for repeat

  return true

end # iszero_numerical

