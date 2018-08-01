%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%    ADI Method for 2D heat equation                                   %
%                                                                      %
%          u_t = u_{xx} + u_{yy} + f(x,t)                              %
%                                                                      %
%    Test problme:                                                     %
%      Exact solution: u(t,x,y) = exp(-t) sin(pi*x) cos(pi*y)          %
%      Source term:  f(t,x,y) = exp(-t) sin(pi*x) cos(pi*y) (2pi^2-1)  %
%                                                                      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 20;
Nt = 20;

h = 1 / Nx;
dt = 1 / Nt;

% uexact = @(x,y,t) exp(-t).*sin(pi*x).*sin(pi*y);
% f = @(x,y,t) exp(-t).*sin(pi*x).*sin(pi*y)*(2*pi*pi-1);

uexact = @(x,y,t) exp(-t).*sin(pi*x).*cos(pi*y);
f = @(x,y,t) exp(-t).*sin(pi*x).*cos(pi*y)*(2*pi*pi-1);

u = zeros(Nx+1,Nx+1,Nt+1);

x = 0:h:1;

xi = h:h:1-h;

t = 0:dt:1;

% Initial condition
u(:,:,1) = uexact(x,x',0);

R = dt/(2*h^2);

e0 = ones(Nx-1,1);
A = spdiags([-R*e0 (1+2*R)*e0 -R*e0], [-1 0 1], Nx-1, Nx-1);
B = spdiags([R*e0 (1-2*R)*e0 R*e0], [-1 0 1], Nx-1, Nx-1);

t0 = 0;
t1 = dt/2;
t2 = dt;
f1 = t1 * f(xi,xi',t1);
f1(1,:) = f1(1,:) + R*uexact(xi,0,t0);
f1(end,:) = f1(end,:) + R*uexact(xi,1,t0);
f2 = t1 * f(xi,xi',t2);
f2(1,:) = f2(1,:) + R*uexact(xi,0,t2);
f2(end,:) = f2(end,:) + R*uexact(xi,1,t2);

u1 =  (B * u(2:end-1,2:end-1,1) + f1) / A;
u(2:end-1,2:end-1,2) = A \ (u1 * B + f2);
u(1,:,2) = uexact(x,0,t2);
u(end,:,2) = uexact(x,1,t2);

u2 = max(max(abs(u(:,:,2) - uexact(x,x',dt))))

figure(1)
surf(u(:,:,1)')
xlabel('x')
ylabel('y')
figure(2)
rotate3d on
surf(u(:,:,2)')
xlabel('x')
ylabel('y')
rotate3d on
