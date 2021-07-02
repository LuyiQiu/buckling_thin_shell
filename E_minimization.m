function S = E_minimization(p,k,a0,alpha,b) % Energy minimization
% use newton's method to calculate energy minimization point a.
% input
% a0 is the initial trial solution
% p,k are pressure and curvature value for energy expression in
% corresponding normalization
% b=1 for balloon and b=0 for shell.

%output
%S.sol is the parameter solution of energy minimization
%S.sod indicates signularity of derivative / det(Hessian matrix)=0.
%S.ni is the number of iteration, when it equals to iteration in the input,
%the number of iteration may not be enough.        
        E = @(a0) E_expression_balloon(p,k,a0,alpha);
        %E = @(a0) E_expression_shell(p,k,a0,alpha);
        S = Newton_method(E,0.00001,50,a0);
end