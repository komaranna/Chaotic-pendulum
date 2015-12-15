function chaotic(k,a0,b0,va0,vb0,m1,m2,l1,l2)
% Chaotic double pendulum in a linearly dissipative field.
% k: linear dissipative constant, negative values result in dissipation.

close all

if (nargin==5)
    %test mode, parameters in SI
    m1=0.2; % masses of pendula
    m2=0.1;
    l1=0.5; % length of pendula
    l2=0.2;
end

g=9.81;
theta1 = m1.*l1.^2./2; % moments of inertia
theta2 = m2.*l2.^2./2;

% max. time[s]
T = 20;
% time step[s]
h=0.001;
N = T./h + 1;

t=linspace(0,T,N); % for plotting the coordinates in time
t2=linspace(h,T,N-1); % for plotting the velocities in time
an=linspace(0,N,N);
bn=linspace(0,N,N);
x1n=linspace(0,N,N);
y1n=linspace(0,N,N);
x2n=linspace(0,N,N);
y2n=linspace(0,N,N);
vx1n=linspace(0,N,N-1);
vy1n=linspace(0,N,N-1);
vx2n=linspace(0,N,N-1);
vy2n=linspace(0,N,N-1);

% initial values
an(1)=a0;
bn(1)=b0;
x1n(1) = l1.*sin(an(1));
y1n(1) = l1.*cos(an(1));
x2n(1) = x1n(1) + l2.*sin(bn(1));
y2n(1) = y1n(1) + l2.*cos(bn(1));
vx1n(1) = va0.*cos(a0);
vy1n(1) = va0.*sin(a0);
vx2n(1) = va0.*cos(a0) + vb0.*cos(b0);
vy2n(1) = va0.*sin(a0) + vb0.*sin(b0);


% first step
an(2) = an(1) + va0.*h;
bn(2) = bn(1) + vb0.*h;
x1n(2) = l1.*sin(an(2));
y1n(2) = l1.*cos(an(2));
x2n(2) = x1n(2) + l2.*sin(bn(2));
y2n(2) = y1n(2) + l2.*cos(bn(2));

% next steps
for i=2:N-1
    an(i+1) = (- k.*an(i-1).*l1.^2./(4.*h) - (theta1+m2.*l1.^2).*an(i-1)./h.^2 + (theta1+m2.*l1.^2).*2.*an(i)./h.^2 - (m1./2+m2).*g.*l2.*sin(an(i))) ./ ((theta1+m2.*l1.^2)./h.^2 - k.*l1.^2./(4.*h));
    bn(i+1) = (- k.*bn(i-1).*l2.^2./(4.*h) - theta2.*bn(i-1)./h.^2 + theta2.*2.*bn(i)./h.^2 - m2./2.*g.*l2.*sin(bn(i))) ./ (theta2./h.^2 - k.*l2.^2./(4.*h));
    
    % calculating coordinates and velocities
    x1n(i+1) = l1.*sin(an(i+1));
    y1n(i+1) = l1.*cos(an(i+1));
    x2n(i+1) = x1n(i+1) + l2.*sin(bn(i+1));
    y2n(i+1) = y1n(i+1) + l2.*cos(bn(i+1));
    
    vx1n(i) = (x1n(i+1) - x1n(i-1))./(2.*h);
    vy1n(i) = (y1n(i+1) - y1n(i-1))./(2.*h);
    vx2n(i) = (x2n(i+1) - x2n(i-1))./(2.*h);
    vy2n(i) = (y2n(i+1) - y2n(i-1))./(2.*h);
end


% Plotting
% Orbits
figure(1)
subplot(121)
plot(x1n,-y1n,'LineWidth',3);
axis square
xlabel('x [m]','FontWeight', 'Bold')
ylabel('y [m]','FontWeight', 'Bold')
title('Orbit of the upper pendulum','FontWeight', 'Bold')

subplot(122)
plot(x2n,-y2n,'LineWidth',3,'Color',[1,0,0]);
axis square
xlabel('x [m]','FontWeight', 'Bold')
ylabel('y [m]','FontWeight', 'Bold')
title('Orbit of the lower pendulum','FontWeight', 'Bold')

% Time evolution of the coordinates
figure(2)
subplot(221)
plot(t,x1n); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('x [m]','FontWeight', 'Bold')
title('Movement of the upper pendulum','FontWeight', 'Bold')

subplot(222)
plot(t,y1n); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('y [m]','FontWeight', 'Bold')

subplot(223)
plot(t2,vx1n); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('v_x [m/s]','FontWeight', 'Bold')

subplot(224)
plot(t2,vy1n); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('v_y [m/s]','FontWeight', 'Bold')


figure(3)
subplot(221)
plot(t,x2n,'Color',[1,0,0]); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('x [m]','FontWeight', 'Bold')
title('Movement of the lower pendulum','FontWeight', 'Bold')

subplot(222)
plot(t,y2n,'Color',[1,0,0]); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('y [m]','FontWeight', 'Bold')

subplot(223)
plot(t2,vx2n,'Color',[1,0,0]); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('v_x [m/s]','FontWeight', 'Bold')

subplot(224)
plot(t2,vy2n,'Color',[1,0,0]); 
xlabel('t [s]','FontWeight', 'Bold')
ylabel('v_y [m/s]','FontWeight', 'Bold')




end




