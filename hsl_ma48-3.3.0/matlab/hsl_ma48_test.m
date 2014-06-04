function hsl_ma48_test()
%
% Unit tests for hsl_ma48 matlab interface
%
fails = 0;

fprintf('Testing toy example:\n')
A = sparse ([1 1 2 2 3 3 3 4 4], [2 4  1 3  1 2 3  1 4], [1.1 3.3, 1.1 4.3, 2.2 4.4 5.5, 3.3 6.6]);
fails = fails + test_with_matrix(A);

fprintf('Testing Wathen(2,4):\n')
A = gallery('wathen', 2, 4);
fails = fails + test_with_matrix(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(fails == 0)
   fprintf('Test OK.\n')
else
   fprintf('Failed %i tests.\n', fails)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fails = test_with_matrix(A)
% Run through all tests with a specified matrix
fails = 0;

control.btf = 1;
if(isreal(A))
   x = rand(size(A,1));
   x2 = rand(size(A,1));
else
   x = rand(size(A,1)) + rand(size(A,1))*i;
   x2 = rand(size(A,1)) + rand(size(A,1))*i;
end
b = A*x;
b2 = A*x2;

fprintf('   - own ordering, separate calls\n')
[handleA, info] = hsl_ma48_factor(A, control);
soln = hsl_ma48_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma48_destroy(handleA);

fprintf('   - colamd ordering, separate calls\n')
P = 1:size(A,1);
Q = colamd(A);
handleA = hsl_ma48_factor(A, control, P, Q);
soln = hsl_ma48_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma48_destroy(handleA);

fprintf('   - colamd ordering, backslash with handle\n')
[soln, info, handleA] = hsl_ma48_backslash(A, b, control, P, Q);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end

fprintf('   - subsequent solve\n')
[soln, info] = hsl_ma48_solve(handleA, b2, control);
res = norm(A*soln - b2, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b2, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma48_destroy(handleA);

fprintf('   - simple backslash\n')
soln = hsl_ma48_backslash(A, b, control);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
