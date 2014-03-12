
% Doc example  Chapter 5.

% Copyright 1999 The MathWorks, Inc.
% $Revision: 1.1.6.1 $

function collect

    y = zeros(1,100); % pre-allocate the matrix 
    for i = 1:100 
        y(i) = collect_one; 
    end

    
    
function y = collect_one
    %#EXTERNAL
    persistent t;
 
    if (isempty(t)) 
        t = 0; 
    else 
        t = t+0.05; 
    end 
    y = sin(t);
