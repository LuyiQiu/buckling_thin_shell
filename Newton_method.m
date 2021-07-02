function S = Newton_method(E,delta,iteration,a0)
% input
% a0 is the initial trial solution
% iteration is the number of iteration 
% delta is the step length to calculate the derivative discretely

%output
%S.sol is the parameter solution of energy minimization
%S.sod indicates signularity of derivative / det(Hessian matrix)=0.
%S.ni is the number of iteration, when it equals to iteration in the input,
%the number of iteration may not be enough.

%the ending criteria : max(step^2) < delta, step= nabla_f/(hessian matrix).

    at=a0; % set the trial solution to a0
    sod=0;
    ni=iteration;
    for iter=1:iteration
        N=[1]; %nabla_f
        for k=1:length(at)
            pr=zeros(1,length(at)); 
            % pr is the direction vector to calculate discrete first order derivative
            % dE/dx_k
            pr(k)=delta; 
            N(k)=(E(at+pr)-E(at))/delta; % N is nabla_f as the discrete derivative
        end
        H=[]; % the hessian matrix
        for k=1:length(at)
            for j=1:length(at)
                pr=zeros(3,length(at));
                % pr is the direction vector to calculate discrete second order derivative
                % d^2E/(dx_kdx_j) or d^2E/(d^2x_k)
                    pr(1,k)=delta;
                    pr(1,j)=delta;
                    pr(2,j)=delta;
                    pr(3,k)=delta;            
                if k==j %d^2E/(d^2x_k)
                    H(k,j)=(E(at+pr(1,:))+E(at-pr(1,:))-2*E(at))/delta^2;
                else  %d^2E/(dx_kdx_j)
                    H(k,j)=(E(at+pr(1,:))-E(at+pr(2,:))-(E(at+pr(3,:))-E(at)))/delta^2;
               end
            end
        end
        if det(H)==0 %the hessian matrix
            sod=1;
            break
        end
        step=(H\N')';
        at=at-step;
        if norm(step) < delta
            ni=iter;
            break
        end
    end
    S.sol=at;
    S.ni=ni;
    S.sod=sod;
end




