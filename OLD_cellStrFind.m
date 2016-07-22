function [index, index_mat] = cellStrFind(str1,str2)

if ~iscell(str1), str1 = num2cell(str1,2); end
if ~iscell(str2), str2 = num2cell(str2,2); end
if size(str1,2) > size(str1,1) str1 = str1'; end
if size(str2,2) > size(str2,1) str2 = str2'; end

str1      = cellfun(@lower,str1,'uniformoutput',false);
str2      = cellfun(@lower,str2,'uniformoutput',false);
index_mat = cell2mat(cellfun(@strcmpi,num2cell(repmat(str1,1,numel(str2)),1),str2','uniformoutput',false));
index     = find(any(index_mat,2));
