function e = eigmatrix(a1)
%EIGMATRIX Returns the eigen value of the given matrix
%    This function returns the eigen value of the input matrix. This
%    function is used to demonstate the functionality of MATLAB Compiler.
%    Refer to the shared library section of MATLAB Compiler for more
%    information.

% Copyright 2003-2007 The MathWorks, Inc.

	try
		%Tries to calculate the eigen value and return it.
		e = eig(a1);
	catch
		%Returns a -1 on error.
		e = -1;
end