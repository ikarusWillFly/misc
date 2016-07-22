function openScripts(path)
if nargin < 1; path = cd; end
files = dir([path,'*.m']);
files = cat(1,{files.name});
cellfun(@(x) edit([path,x]),files);
end