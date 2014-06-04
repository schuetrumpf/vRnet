function hsl_ma48_expert(handle)
% HSL_MA48_DESTROY  Free memory associated with factorization.
%     hsl_ma48_destroy(handle) will free all memory associated with handle.
%     The factorization may not be reused again.
%
%     Usage: hsl_ma48_destroy(handle)
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] MA48, a Fortran code for direct colution of sparse unsymmetric linear
%         systems of equations. I.S. Duff and J.K. Reid. Report RAL-93-072.
%
%     See also: ma48_backslash, ma48_factor, ma48_solve

hsl_ma48_expert('destroy', handle)
