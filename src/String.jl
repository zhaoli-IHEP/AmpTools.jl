#############################################################
#
# Although they do not have relation with SymEngine,
#   the following functions have been widely used in our projects.
#
#############################################################


#####################################################################
"""
    box_message( message::String; color=:light_cyan )::Nothing

Return nothing.
"""
function box_message( message::String; color=:light_cyan )::Nothing
#####################################################################

  line_list = filter( !isempty, split( message, "\n" ) )

  len = (maximum∘map)( length, line_list ) 

  line_list = map( x -> "[ $(x)$(" "^(len-length(x))) ]", line_list )

  bar = "-"^(len+4)

  new_message = join( line_list, "\n" )

  printstyled( """

  $bar
  $(new_message)
  $bar

  """, color=color )

  return nothing

end # function box_message






######################################
function seq_replace(
    str::String,
    rule_dict::Dict{String,String}
)::String
######################################

  for one_rule in collect( rule_dict )
    str = replace( str, one_rule )
  end # for one_rule

  return str

end # function seq_replace 

######################################
function seq_replace(
    str::String,
    rule_list::Vector{Pair{String,String}}
)::String
######################################

  for one_rule in rule_list 
    str = replace( str, one_rule )
  end # for one_rule

  return str

end # function seq_replace 

######################################
function seq_replace(
    str::String,
    rule_list...
)::String
######################################

  for one_rule in rule_list 
    str = replace( str, one_rule )
  end # for one_rule

  return str

end # function seq_replace 




############################################
@inline add_quote( str::Any )::String = "\"$str\""
############################################
@inline add_quote( str_list::Vector{Any} )::Vector{String} = add_quote.(str_list)
############################################


###############################################
# backup before mkdir
function bk_mkdir( dir_name::String )::Nothing
###############################################
  
  if isdir( dir_name )
    mv( dir_name, "$(dir_name)_$(now())" )
  end # if
  mkdir( dir_name )

  return nothing

end # function bk_mkdir




