function r = realacos(x)
% REALACOS Arc-cosine, output guaranteed to be non-complex
%
% This function generates an error if the input has absolute value greater 
% than 1.0, because only those values between -1.0 and 1.0 have non-complex
% arc cosine.

% Copyright 2005 The MathWorks, Inc.

    if x > 1.0
        error(['Input X must be <= 1.0. (Was: ' num2str(x) ')']);
    end

    r = acos(x);

