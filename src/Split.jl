
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



