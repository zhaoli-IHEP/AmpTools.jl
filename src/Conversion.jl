
###########################################
gen_sorted_str( ::Val{:Symbol}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Integer}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Rational}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Complex}, expr::Basic )::String = string(expr)
###########################################


###########################################
function gen_sorted_str( ::Val{:Pow}, expr::Basic )::String
###########################################
  arglist = get_args(expr)
  return ":Pow( $(gen_sorted_str(arglist[1])), $(arglist[2]) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:FunctionSymbol}, expr::Basic )::String
###########################################
  name = replace( get_name(expr), "Trace" => "DiracTrace" )
  return "$(name)( $(join( map( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:Mul}, expr::Basic )::String
###########################################
  return ":Mul( $(join( (sort∘map)( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:Add}, expr::Basic )::String
###########################################
  return ":Add( $(join( (sort∘map)( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str
 
 
###########################################
"""
    gen_sorted_str( expr::Basic )::String
 
This is a generic interface for different classes of the `expr`.
And it will generate the sorted string format for the expression `expr`.
"""
function gen_sorted_str( expr::Basic )::String
###########################################
 
  expr_class = SymEngine.get_symengine_class(expr)
  return gen_sorted_str( Val(expr_class), expr )
 
end # function gen_sorted_str
 


###############################
function to_String_dict(
    dict::Dict{Basic,Basic}
)::Dict{String,String}
###############################

  if isempty(dict)
    return Dict{String,String}()
  end # if

  return (Dict∘map)( x -> string(x[1]) => string(x[2]), collect(dict) )

end # function to_String_dict

###############################
function to_Basic_dict(
    dict::Dict{String,String}
)::Dict{Basic,Basic}
###############################

  if isempty(dict)
    return Dict{Basic,Basic}()
  end # if

  return (Dict∘map)( x -> Basic(x[1]) => Basic(x[2]), collect(dict) )

end # function to_Basic_dict


###############################
function to_Basic_dict(
    dict::Dict{Any,Any}
)::Dict{Basic,Basic}
###############################

  return (Dict∘map)( x -> (Basic∘string)(x[1]) => (Basic∘string)(x[2]), collect(dict) )

end # function to_Basic_dict



################################
function to_Basic_dict(
    dict_list::Union{ Vector{Vector{String}}, Vector{Tuple{String,String}} }
)::Dict{Basic,Basic}
################################

  result_dict = Dict{Basic,Basic}()
  for one_pair in dict_list
    @assert length(one_pair) == 2
    push!( result_dict, Basic(one_pair[1]) => Basic(one_pair[2]) )
  end # for one_pair

  return result_dict

end # function to_Basic_dict



##############################################################
@inline to_Basic( str_list::Vector{String} )::Vector{Basic} = map( Basic, str_list )
##############################################
@inline to_String( ex_list::Vector{Basic} )::Vector{String} = map( string, ex_list )
##############################################


###########################################
gen_mma_str( ::Val{:Symbol}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Integer}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Rational}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Complex}, expr::Basic )::String = string(expr)
###########################################


###########################################
function gen_mma_str( ::Val{:Pow}, expr::Basic )::String
###########################################
  arglist = get_args(expr)
  return "Power[ $(gen_mma_str(arglist[1])), $(arglist[2]) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:FunctionSymbol}, expr::Basic )::String
###########################################
  name = replace( get_name(expr), "Trace" => "DiracTrace" )
  return "$(name)[ $(join( map( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:Mul}, expr::Basic )::String
###########################################
  return "Times[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:Add}, expr::Basic )::String
###########################################
  return "Plus[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str


###########################################
"""
    gen_mma_str( expr::Basic )::String

This is a generic interface for different classes of the `expr`.
And it will generate the mathematica code for the expression `expr`.
"""
function gen_mma_str( expr::Basic )::String
###########################################

  expr_class = SymEngine.get_symengine_class(expr)
  return gen_mma_str( Val(expr_class), expr )

end # function gen_mma_str




