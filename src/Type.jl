



####################
function is_FunctionSymbol(
    expr::Basic
)::Bool
####################

  return SymEngine.get_symengine_class(expr) == :FunctionSymbol

end # function is_FunctionSymbol


