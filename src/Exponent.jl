

###################################################################
function get_exponent( ::Val{:Pow}, poly::Basic, x::Basic, op::Symbol = :Max )::Int64
###################################################################
  arg_list = get_args(poly)
  #return arg_list[1] == x ? Int64( arg_list[2] ) : 0
  return get_exponent( arg_list[1], x, op )*Int64( arg_list[2] )

end # function get_exponent


###################################################################
function get_exponent( ::Val{:Symbol}, poly::Basic, x::Basic, op::Symbol = :Max )::Int64
###################################################################
  return poly == x ? 1 : 0
end # function get_exponent

get_exponent( ::Val{:FunctionSymbol}, poly::Basic, x::Basic, op::Symbol = :Max)::Int64 = 0
get_exponent( ::Val{:Integer}, poly::Basic, x::Basic, op::Symbol = :Max)::Int64 = 0
get_exponent( ::Val{:Rational}, poly::Basic, x::Basic, op::Symbol = :Max)::Int64 = 0
get_exponent( ::Val{:RealMPFR}, poly::Basic, x::Basic, op::Symbol = :Max)::Int64 = 0


###################################################################
function get_exponent( ::Val{:Complex}, poly::Basic, x::Basic, op::Symbol = :Max)::Int64
###################################################################
  return x == Basic(im) ? 1 : 0
end # function get_exponent


###################################################################
function get_exponent( ::Val{:Add}, poly::Basic, x::Basic, op::Symbol = :Max )::Int64
###################################################################
  if op == :Max
    return max( map( p_ -> get_exponent(p_,x,op), get_args(poly) )... )
  elseif op == :Min
    return min( map( p_ -> get_exponent(p_,x,op), get_args(poly) )... )
  else
    error( "Exception" )
  end # if
end # function get_exponent

###################################################################
function get_exponent( ::Val{:Mul}, poly::Basic, x::Basic, op::Symbol = :Max )::Int64
###################################################################
  return (sum∘map)( p_ -> get_exponent(p_,x,op), get_args(poly) )
end # function get_exponent



###################################################################
"""
    get_exponent( poly::Basic, x::Basic, op::Symbol = :Max )::Int64

Extract the degree (the max exponent by default) of `x` in polynomial `poly`.
"""
function get_exponent( poly::Basic, x::Basic, op::Symbol = :Max )::Int64
###################################################################

  #expanded_poly = expand(poly)
  #poly_class = SymEngine.get_symengine_class( expanded_poly )
  #return get_exponent( Val(poly_class), expanded_poly, x, op )

  if iszero(poly)
    if op == :Max
      return -99999999
    elseif op == :Min
      return 99999999
    else
      error("Exception")
    end # if
  end # if

  poly_class = SymEngine.get_symengine_class( poly )
  return get_exponent( Val(poly_class), poly, x, op )

end # function get_exponent



###################################################################
"""
    get_degree( poly::Basic, x::Basic )::Int64

Extract the degree (the max exponent) of `x` in polynomial `poly`.
"""
function get_degree( poly::Basic, x::Basic )::Int64
###################################################################

 #expanded_poly = expand(poly)
 #poly_class = SymEngine.get_symengine_class( expanded_poly )
 #return get_exponent( Val(poly_class), expanded_poly, x, :Max )

  poly_class = SymEngine.get_symengine_class( poly )
  return get_exponent( Val(poly_class), poly, x, :Max )

end # function get_degree

