%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function output = find_outlier( input, alpha, k)
% Rosner, Bernard. „Percentage Points for a Generalized ESD Many-Outlier Procedure“. Technometrics 25, Nr. 2 (1. Mai 1983): 165–72.
%
% written by R.Bauer for CIN AG NPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = find_outlier( input, alpha, k)

if size(input,1) > size(input,2), input = input'; end % turn row to column

posnan      = find(isnan(input));
pos         = [];
N           = sum(~isnan(input));
if nargin < 2, alpha = 0.05; end
if nargin < 3, k = floor(N/10); end 

R = zeros(k+1, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculations                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

M = nanmean(input); 
[value, index] = sort(abs(input - M));

for i = 0:k, 
    V = value(1:N-i); 
    R(i+1) = abs(V(N-i)-nanmean(V)) / nanstd(V); 
end; 

for i = 1:k,  
    crit        = 1-(alpha/((2*(k-i+1))));
    t           = tinv(crit, N-i-1);
    lambda(i)   =(N-i)*t./sqrt(((N-i-1+t^2)*(N-i)));
    if R(i)     > lambda (i)
        pos  = index(N-i+1:end);   
    end
    
    cleared = input;
    
    if ~isempty(pos),
        cleared(pos) = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% output                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
output.excl     = sort(unique([posnan,pos]));
output.incl     = 1:length(input); output.incl(pos) = [];
output.cleared  = cleared;
end
