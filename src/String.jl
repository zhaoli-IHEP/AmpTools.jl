
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
function convert_to_String_dict(
    dict::Dict{Basic,Basic}
)::Dict{String,String}
###############################

  return (Dict∘map)( x -> string(x[1]) => string(x[2]), collect(dict) )

end # function convert_to_String_dict

###############################
function convert_to_Basic_dict(
    dict::Dict{String,String}
)::Dict{Basic,Basic}
###############################

  return (Dict∘map)( x -> Basic(x[1]) => Basic(x[2]), collect(dict) )

end # function convert_to_Basic_dict

