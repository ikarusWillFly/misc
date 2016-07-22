 files = dir('*.m');
 files = cat(1,{files.name});
 cellfun(@edit,files)