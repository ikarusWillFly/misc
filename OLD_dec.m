function rounded = dec(num,decimal)
if nargin < 2, decimal = 2; end
factor = (10^decimal);
rounded = round(num*factor)/factor;
end

