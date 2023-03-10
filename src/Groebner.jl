

######################################
# https://math.colorado.edu/~rohi1040/expository/groebner_bases.pdf
######################################


################################
function get_monomial_weight( 
    monomial::Basic, 
    xi_list::Vector{Basic} 
)::Vector{Int64}
################################

  the_class = SymEngine.get_symengine_class(monomial)
  if the_class in [:Integer,:Rational,:Complex]
    return zeros(Int64,length(xi_list)+1)
  end # if

  @assert the_class in [:Mul,:Pow,:Symbol]
  @assert free_symbols(monomial) ⊆ xi_list

  xpt_list = map( xi -> get_exponent(monomial,xi), xi_list )
  # graded lexicographic order
  #return vcat( [sum(xpt_list)], xpt_list )
  return vcat( [zero(Basic)], xpt_list )

end # function get_monomial_weight

#############################
# leading term
function get_LT( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return polynomial_expand
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  return term_list[pos]

end # function get_LT


#############################
# leading monomial e.g. power product
function get_LM( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return polynomial_expand
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  xpt_list = weight[2:end]
  return prod( xi_list .^ xpt_list )

end # function get_LM

#############################
# leading coefficient 
function get_LC( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return polynomial_expand
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  return subs( term_list[pos], xi_list .=> one(Basic) )

end # function get_LC


####################################################
# Define the function to compute the S-polynomial
function get_spoly( 
    f::Basic, 
    g::Basic, 
    xi_list::Vector{Basic} 
)::Basic
####################################################

  f_LT = get_LT( f, xi_list )
  g_LT = get_LT( g, xi_list )

  f_LT_xpt_list = get_monomial_weight( f_LT, xi_list )[2:end]
  g_LT_xpt_list = get_monomial_weight( g_LT, xi_list )[2:end]

  lcm_xpt_list = max.( f_LT_xpt_list, g_LT_xpt_list )
  lcm_fg = prod( xi_list .^ lcm_xpt_list )

  return expand( lcm_fg / f_LT * f - lcm_fg / g_LT * g )

end # function get_spoly



###########################################################
# Define the function to reduce a polynomial with respect to a Gröbner basis
function reduce_Groebner( 
    f::Basic, 
    G_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Basic
###########################################################

  r = copy(f)
  for G in G_list

    LT_g = get_LT( G, xi_list )
    weight_LT_g = get_monomial_weight( LT_g, xi_list )
    while true

      LT_r = get_LT( r, xi_list )
      weight_LT_r = get_monomial_weight( LT_r, xi_list )

      if all( weight_LT_r .>= weight_LT_g ) 
        r = expand( r - LT_r / LT_g * G )
        if iszero(r)
          return r
        end # if
      else
        break
      end # if

    end # while

  end # for G

  return r

end # function reduce_Groebner





###########################################################
# Define the function to reduce a polynomial with respect to a Gröbner basis
function reduce_Groebner_v2( 
    f::Basic, 
    G_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Basic
###########################################################

  r = copy(f)
  while true
    old_r = copy(r)

    for G in G_list

      LT_g = get_LT( G, xi_list )
      weight_LT_g = get_monomial_weight( LT_g, xi_list )

      while true

        LT_r = get_LT( r, xi_list )
        weight_LT_r = get_monomial_weight( LT_r, xi_list )
  
        if all( weight_LT_r .>= weight_LT_g ) 
          r = expand( r - LT_r / LT_g * G )
          if iszero(r)
            return zero(Basic)
          end # if
        else
          break
        end # if

      end # while

    end # for G

    if r == old_r 
      break
    end # if

  end # while

  return r

end # function reduce_Groebner_v2



















################################################################################
# Define the function to compute the Gröbner basis using the Buchberger's algorithm
function get_Groebner_basis( 
    F_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Vector{Basic}
################################################################################

  G_list = copy(F_list)

  while true
    Gp_list = copy(G_list)
    for (p,q) in (collect∘powerset)( Gp_list, 2, 2 )
      s = get_spoly( p, q, xi_list)
      s = reduce_Groebner_v2( s, Gp_list, xi_list )
      if !iszero(s)
#@show p q s G_list
        push!( G_list, s )
      end # if
    end # for one_pair

    if length(G_list) == length(Gp_list) &&
       (isempty∘setdiff)( G_list, Gp_list )
      break
    end # if
  end # while 

  return G_list

end # function get_Groebner_basis




######################################
# Interface to Groebner.jl
function get_Groebner_basis_v2(
    F_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Vector{Basic}
######################################

  xi_str_list = map( string, xi_list )
  R, new_xi_list = PolynomialRing(QQ, xi_str_list);

  polys = Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}()
  for F in F_list 
    term_list = get_add_vector_expand(F)
    new_F = zero(R)
    for one_term in term_list
      the_coeff = convert( Int64, subs( one_term, Dict(xi_list .=> one(Basic)) ) )
      xpt_list = get_monomial_weight( one_term, xi_list )[2:end]
      new_F += the_coeff * prod( new_xi_list .^ xpt_list )
    end # for one_term
    push!( polys, new_F )
  end # for F

  G_list = groebner(polys)

  new_G_list = Vector{Basic}()
  for G in G_list
    new_G = zero(Basic)
    n_term = length(G)
    for index in 1:n_term
      the_coeff = AbstractAlgebra.coeff(G,index) 
      one_term = term(G,index)
      xpt_list = map( x->AbstractAlgebra.degree(one_term,x), new_xi_list )

      new_term = the_coeff * prod( xi_list .^ xpt_list )
      new_G += new_term
    end # for index
    push!( new_G_list, new_G )
  end # for G

  return new_G_list

end # function get_Groebner_basis_v2

















