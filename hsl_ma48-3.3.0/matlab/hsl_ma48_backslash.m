function [X, varargout] = hsl_ma48_backslash(A, B, varargin)
% HSL_MA48_BACKSLASH Sparse Symmetric Indefinite Linear Solve.
%     X = hsl_ma48_backslash(A, B) solves the equation AX=B for X by means of
%     an unsymmetric factorization PAQ = LU. P is generated
%     automatically to reduce fill-in.
%
%     Usage: X = hsl_ma48_backslash(A, B)
%            [X, info, handle] = ma48_backslash(A, B, control, P)
%
%     control is a structure described below. P is a permuaton such as that
%     output from symamd(A). info is a structure described below. handle
%     provides a reference to the factorization so it may be reused at a later
%     date.
% 
%     control may have the following components set. If they are not set then
%     the stated default is used.
%     control.hermitian    - True or false. Determines if a complex matrix is
%                            treated as Hermitian (true) or symmetric (false).
%                            Default is false.
%     control.nb           - Block size to be used. Default is 256.
%     control.nemin        - Maximum number of columns in candidates for
%                            supernode amalgamation. Default is 32.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%     control.scaling      - Determines if scaling is to be used with values:
%                                   1 : MC77 in the one norm
%                            otherwise:  no scaling
%                            Default is 1.
%     control.small        - Pivots of modulus less than this are treated as
%                            zero. Default is 1e-20.
%     control.static       - If greater than zero static pivoting is used.
%                            Default is 0.0.
%     control.u            - Initial relative pivot tolerance threshold. Default
%                            is 0.01.
%     control.umin         - Relaxed relative pivot tolerance threshold. Default
%                            is 1.0.
%
%     On return, info will have the following components set.
%     info.matrix_rank        - Number of non-zero pivots.
%     info.num_delay          - Number of delayed pivots.
%     info.num_factor         - Number of entries in the factors (after
%                               supernode amalgamation and pivoting).
%     info.num_flops          - Number of floating point operations to form
%                               factors (after supernode amalgamation and
%                               pivoting).
%     info.num_neg            - Number of negative pivots in factors.
%     info.num_perturbed      - Number of perturbed pivots when static pivoting
%                               is used.
%     info.num_two            - Number of 2x2 pivots used in factorization.
%     info.order              - Ordering used. One of 'MeTiS', 'AMD' or 'user'.
%                               If a user has provided a permutation that is
%                               used, otherwise a MeTiS ordering is used.
%                               Should MeTiS not be available then AMD is used
%                               instead.
%     info.usmall             - Threshold parameter actually used. This will
%                               only differ from control.u if control.umin <
%                               control.u.
%     info.order_time         - Wall clock time for Fortran ordering routine
%     info.analyse_time       - Wall clock time for Fortran ma48_analyse call
%     info.factor_solve_time  - Wall clock time for Fortran ma48_factor_solve
%                               call
%
%     Note that if the optional return value handle is used, then memory will
%     be tied up storing the factorization. This may be recovered by calling
%     ma48_destroy.
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] MA48, a Fortran code for direct colution of sparse unsymmetric linear
%         systems of equations. I.S. Duff and J.K. Reid. Report RAL-93-072.
%
%     See also: ma48_destroy, ma48_factor, ma48_solve

% Note: ma48 doesn't offer a combined factor_solve operation, so we just do
% a factor then a solve before destroying the intermediate data (unless handle
% is requested to be returned)

optargin = size(varargin,2);
optargout = max(nargout,1)-1;
if(optargin>3)
   error ('Too many arguments')
end
if(optargout == 0)
   if(optargin == 0)
      handle = hsl_ma48_expert('factor', A);
      X = hsl_ma48_expert('solve', handle, B);
   elseif(optargin == 1)
      handle = hsl_ma48_expert('factor', A, varargin{1});
      X = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 2)
      handle = hsl_ma48_expert('factor', A, varargin{1}, varargin{2});
      X = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 3)
      handle = ...
         hsl_ma48_expert('factor', A, varargin{1}, varargin{2}, varargin{3});
      X = hsl_ma48_expert('solve', handle, B, varargin{1});
   end
   hsl_ma48_expert('destroy', handle);
elseif(optargout == 1)
   if(optargin == 0)
      [handle, info] = hsl_ma48_expert('factor', A);
      [X, infoS] = hsl_ma48_expert('solve', handle, B);
   elseif(optargin == 1)
      [handle, info] = hsl_ma48_expert('factor', A, varargin{1});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 2)
      [handle, info] = hsl_ma48_expert('factor', A, varargin{1}, varargin{2});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 3)
      [handle, info] = ...
         hsl_ma48_expert('factor', A, varargin{1}, varargin{2}, varargin{3});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   end
   hsl_ma48_expert('destroy', handle);
   info.solve_time = infoS.solve_time;
   varargout(1) = {info};
elseif(optargout == 2)
   if(optargin == 0)
      [handle, info] = hsl_ma48_expert('factor', A);
      [X, infoS] = hsl_ma48_expert('solve', handle, B);
   elseif(optargin == 1)
      [handle, info] = hsl_ma48_expert('factor', A, varargin{1});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 2)
      [handle, info] = hsl_ma48_expert('factor', A, varargin{1}, varargin{2});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   elseif(optargin == 3)
      [handle, info] = ...
         hsl_ma48_expert('factor', A, varargin{1}, varargin{2}, varargin{3});
      [X, infoS] = hsl_ma48_expert('solve', handle, B, varargin{1});
   end
   info.solve_time = infoS.solve_time;
   varargout(1) = {info};
   varargout(2) = {handle};
else
    error ('Too many output arguments')
end
