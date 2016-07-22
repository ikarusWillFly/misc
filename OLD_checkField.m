function structure = checkField(structure,field,value)
% insert the value to a non existent field 
if ~isfield(structure,field), structure.(field) = value; end