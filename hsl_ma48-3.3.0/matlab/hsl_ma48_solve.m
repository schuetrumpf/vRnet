function [X, info] = hsl_ma48_solve(handle, B, varargin)
% HSL_MA48_SOLVE  Sparse Symmetric Indefinite Solve.
%     X = hsl_ma48_solve(handle, B) solves the equation AX=B for X given
%     precomputed factors associated with handle. The handle must have been
%     obtained by a prior call to hsl_ma48_factor of hsl_ma48_backslash.
%
%     Usage: X = hsl_ma48_solve(handle, B)
%            [X, info] = hsl_ma48_solve(handle, B, control)
%
%     The optional argument CONTROL may have the following components set. If
%     they are not set then the stated default is used.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%
%     The optional return value INFO will have some of the following components
%     set on exit.
%     info.solve_time         - Wall clock time for Fortran ma48_solve call
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] MA48, a Fortran code for direct colution of sparse unsymmetric linear
%         systems of equations. I.S. Duff and J.K. Reid. Report RAL-93-072.
%
%     See also: ma48_backslash, ma48_destroy, ma48_factor

optargin = size(varargin,2);
if(optargin == 0)
   [X, info] = hsl_ma48_expert('solve', handle, B);
elseif(optargin == 1)
   [X, info] = hsl_ma48_expert('solve', handle, B, varargin{1});
else
   error ('Too many arguments')
end
