% http://blogs.mathworks.com/community/2012/02/13/parsing-inputs/
function inputs = parseinputtest(varargin)
p = inputParser;
p.addRequired('val',@isscalar);
%p.addOptional('ntimes',1,@isscalar);
p.addParamValue('title','Default title',@isstr);

p.parse(varargin{:});
 inputs = p.Results;
 
 % Invoke: parseinputtest(1,'title','test')