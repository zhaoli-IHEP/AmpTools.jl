
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





#########################
function findmatched(
    str::Union{ String, SubString{String} }, 
    start_pos::Int64
)::Int64
#########################

  left = str[start_pos]
  if left == '('
    right = ')' 
    rx = r"[()]"
  elseif left == '['
    right = ']' 
    rx = r"[\[\]]"
  elseif left == '{'
    right = '}' 
    rx = r"[{}]"
  else
    error("Exception")
  end # if

  pos = (first∘findnext)( rx, str, start_pos )
  lr = str[pos]
  
  count = 1
  while true
    pos = (first∘findnext)( rx, str, pos+1 )
    lr = str[pos]

    if lr == left
      count += 1
    elseif lr == right
      count -= 1
    else
      error("Exception")
    end # if

    if iszero(count)
      return pos
    end # if
  end # while

end # function findmatched


