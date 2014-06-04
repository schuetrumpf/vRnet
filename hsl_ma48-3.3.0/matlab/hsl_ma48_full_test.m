function hsl_ma48_full_test

clear all;

A = gallery('poisson', 2); % yes its spd, not unsymmetric - shouldn't matter
x = [1; 2; 3; 4];
b = A*x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insufficient arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
   handle = hsl_ma48_factor();
   error('Unexpected success at insufficient inputs to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "A" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_factor')
   end
end

try
   X = hsl_ma48_backslash();
   error('Unexpected success at insufficient inputs to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "A" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_backslash')
   end
end

try
   X = hsl_ma48_backslash(A);
   error('Unexpected success at insufficient inputs to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "B" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_backslash')
   end
end

try
   X = hsl_ma48_solve();
   error('Unexpected success at insufficient inputs to hsl_ma48_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "handle" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_solve')
   end
end

handle = hsl_ma48_factor(A);
try
   X = hsl_ma48_solve(handle);
   error('Unexpected success at insufficient inputs to hsl_ma48_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "B" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_solve')
   end
end
hsl_ma48_destroy(handle);

try
   hsl_ma48_destroy();
   error('Unexpected success at insufficient inputs to hsl_ma48_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "handle" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma48_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Too many input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

control.btf=1;
P = 1:size(A,1);
Q = colamd(A);
extra = 1;

try
   handle = hsl_ma48_factor(A, control, P, Q, extra);
   error('Unexpected success at too many inputs to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma48_factor')
   end
end

handle = hsl_ma48_factor(A);
try
   X = hsl_ma48_solve(handle, b, control, extra);
   error('Unexpected success at too many inputs to hsl_ma48_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma48_solve')
   end
end
hsl_ma48_destroy(handle);

try
   X = hsl_ma48_backslash(A, b, control, P, Q, extra);
   error('Unexpected success at too many inputs to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma48_backslash')
   end
end

try
   hsl_ma48_destroy(handle, extra);
   error('Unexpected success at too many inputs to hsl_ma48_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many input arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma48_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Too many output arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
   [handle, info, extra] = hsl_ma48_factor(A, control, P);
   error('Unexpected success at too many outputs to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma48_factor')
   end
end

handle = hsl_ma48_factor(A);
try
   [x, info, extra] = hsl_ma48_solve(handle, b);
   error('Unexpected success at too many outputs to hsl_ma48_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma48_solve')
   end
end
hsl_ma48_destroy(handle);

try
   [x, info, handle, extra] = hsl_ma48_backslash(A, b);
   error('Unexpected success at too many outputs to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma48_backslash')
   end
end

try
   extra = hsl_ma48_destroy(handle);
   error('Unexpected success at too many outputs to hsl_ma48_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma48_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not a sparse matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adense = full(A);

try
   [handle, info] = hsl_ma48_factor(Adense);
   error('Unexpected success at dense matrix to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Error in argument A. Expected sparse matrix.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at dense matrix to hsl_ma48_factor')
   end
end

try
   [handle, info] = hsl_ma48_backslash(Adense, b);
   error('Unexpected success at dense matrix to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Error in argument A. Expected sparse matrix.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at dense matrix to hsl_ma48_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A not square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Arect = [A [1; 2; 3; 4;]];

try
   [handle, info] = hsl_ma48_factor(Arect);
   error('Unexpected success at rectangular matrix to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('The matrix must be square');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at rectangular matrix to hsl_ma48_factor')
   end
end

try
   x = hsl_ma48_backslash(Arect, b);
   error('Unexpected success at rectangular matrix to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('The matrix must be square');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at rectangular matrix to hsl_ma48_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A and b inconsistent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b5 = [1; 2; 3; 4; 5];

handle = hsl_ma48_factor(A);
try
   soln = hsl_ma48_solve(handle, b5);
   error('Unexpected success at inconsistent A, b to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Dimensions of A and b inconsistent: A=    4x    4, b=    5x    1');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at inconsistent A, b to hsl_ma48_factor')
   end
end
hsl_ma48_destroy(handle);

try
   soln = hsl_ma48_backslash(A, b5);
   error('Unexpected success at inconsistent A, b to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Dimensions of A and b inconsistent: A=    4x    4, b=    5x    1');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at inconsistent A, b to hsl_ma48_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pbad = [1 2 2 3];

try
   handle = hsl_ma48_factor(A, control, Pbad, Q);
   error('Unexpected success at bad P to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Problem with P or Q.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma48_factor')
   end
end

try
   x = hsl_ma48_backslash(A, b, control, Pbad, Q);
   error('Unexpected success at bad P to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Problem with P or Q.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma48_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q absent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
   handle = hsl_ma48_factor(A, control, P);
   error('Unexpected success at bad P to hsl_ma48_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('If P is present, Q must also be present.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma48_factor')
   end
end

try
   x = hsl_ma48_backslash(A, b, control, P);
   error('Unexpected success at bad P to hsl_ma48_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('If P is present, Q must also be present.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma48_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('All tests succeeded.\n')
