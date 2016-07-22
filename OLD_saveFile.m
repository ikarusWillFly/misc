function fileName = saveFile(data, path, fileName, extension)
if nargin < 1, data      = [];     end
if nargin < 2, path      = cd;     end
if nargin < 3, fileName  = 'test'; end
if nargin < 4, extension = 'mat';  end
% saves the data with unique file names
nFile    = 1;
filename = @(task,nFile) [task,'_',num2str(nFile),'.',extension];

files = dir([path,fileName,'*']);
files = {files.name};
while any(strcmpi(files,filename(fileName,nFile)))
    nFile = nFile + 1;
end
fileName = [path,fileName,'_',num2str(nFile)];
if ~isempty(data)
save(fileName,'data')
disp(['file succesfully saved in ',fileName])
end

