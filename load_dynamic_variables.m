function [varNameOut,message] = load_dynamic_variables(data,variable_name_out,variable_name_in)
% [varNameOut,message] = load_dynamic_variables(varNameIn,varNameOut)
% load and name variables dynamically
%
% if DATA is a string, the function loads the variable VARIABLE_NAME_IN
% from the path DATA and renames it VARIABLE_NAME_OUT
%
% else, if DATA is the DATA itself, it gets renamed VARIABLE_NAME_OUT
%
% the output MESSAGE contains a warning of the action performed by the
% function
if ischar(data);
    tmp   = load(data,variable_name_in);
    eval([variable_name_out,' = tmp.(variable_name_in);']);
    message = sprintf('LOADING\t%s\tin\t%s',variable_name_out,variable_name_in);  
else
    eval([variable_name_out,' = variable_name_in;'])
    message = sprintf('USING INPUT %s',variable_name_out);
end
end