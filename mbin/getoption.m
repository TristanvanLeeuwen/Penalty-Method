function s = getoption(options,field,default)
% Parse option from struct
%
% use:
%   s = getoption(options,field,default)
%
% input:
%   options - struct
%   field   - name options to be parsed
%   default - default value
%
% output:
%   s - value of option

s = default;
if isfield(options,field)
    s = getfield(options,field);
end
