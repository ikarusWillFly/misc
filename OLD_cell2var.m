function cell = cell2var(var)

if iscell(var)
    cell = cat(1,var{:});
else 
    cell = var;
end