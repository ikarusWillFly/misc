% output = smooth2D(input,window,mode,padding)
% window should be specified as [dim1,dim2]
% mode 1: boxcar
% mode 2: gauss
% mode 3: parzen
% mode 4: hanning

function output = smooth2D(input,window,mode,padding)
    
    if nargin <3, mode = 1; end
    if nargin <4, padding = 'replicate'; end
    if mode == 1, 
        FILT = ones(window(1),1)*ones(window(2),1)';        
        FILT = FILT./numel(FILT);
    elseif mode == 2, 
       FILT =  gausswin(window(1))*gausswin(window(2))';
    elseif mode == 3, 
       FILT =  parzenwin(window(1))*parzenwin(window(2))';
    elseif mode == 4,
       FILT =  hanning(window(1))*hanning(window(2))';    
    elseif mode == 5,
       FILT =  boxcar(window(1))*parzenwin(window(2))'./window(1);
    end
    FILT    = FILT./sum(sum(FILT));
    if ~strcmp(padding,'none'),
        output  = convn(padarray(input,window./2,padding),FILT,'same');
    else
        output  = convn(input,FILT,'same');
    end
    output  = output(1+(window(1)./2):end-window(1)/2,1+(window(2)./2):end-window(2)/2);

end