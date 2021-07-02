function Energy = E_expression_balloon (p,k,a,alpha) % Energy expression 

    Energy=0; % initiation
    phi=[];
    x=[];
    y=[];
    e0=[];
    
    dt=0.001;
    theta=0:dt:2*pi; 
    
    % phi=a(1)*sin(t)+a(2)*sin(2*t)+a(3)*sin(3*t)+a(4)*sin(4*t); %rotation angle function phi
    syms t
    phi = @(t,a) a(1).*sin(t)+a(2).*sin(2.*t)+a(3).*sin(3.*t)+a(4).*sin(4.*t);
    dphi = @(t,a) a(1).*cos(t)+a(2).*2.*cos(2*t)+a(3).*3.*cos(3*t)+a(4).*4.*cos(4*t);
    
    for i=1:length(theta)
        x(i) = integral(@(t) cos(t+phi(t,a)), 0, dt*i);
        y(i) = integral(@(t) sin(t+phi(t,a)), 0, dt*i);
    end

    D_area=sum(x.*sin(theta+phi(theta,a)))*dt-pi; %dimensionless parameter A_bar
    e0=(p*D_area)/(2*pi) - k/(2*pi)*sum(y)*dt; %dimensionless parameter e0_bar

    Energy =0.5*alpha^2*integral (@(t) dphi(t,a).*dphi(t,a), 0, 2*pi) ...
    + 0.5*sum((e0+k.*y).^2)*dt - p * sum((1+e0+k*y).*x.*sin(theta+phi(theta,a)))*dt ...
    + p * pi * (1+e0+k*sum(y)*dt/(2*pi))...
    + a(5)*x(end);
end


