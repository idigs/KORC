function [x,y,z] = torusmap(r, R, v)
% xyz = torusmap(r, R, v)
% INPUTS:
% R large radius, r, cross section radius of the torus
% v is 3 x n array represent n vectors of uniform random in the cube (0,1)^3.
% OUTPUT:
% xyz is 3 x n array 3D coordinates mapped uniformly in the torus
a = r/R;
t = (2/3) * a;
lo = -a-2;
hi = -a+2;

% Intialize by histogram
m = 9;
knots = linspace(0,pi,m);
c = cos(knots);
c2 = c.^2;
s2 = 1 - c2;
s = sqrt(s2);
C2 = (t*s2-c).*s + knots;
u = pi * v(1,:);
phi = interp1(C2, knots, u, 'linear', 'extrap');
c = cos(phi);
q = -(a*c + 2).*c;
q = max(min(q,hi),lo);

niter = 10;
for k = 1:niter
    c = (sqrt(1- a*q)-1)/a;
    c2 = c.^2;
    s2 = 1 - c2;
    s = sqrt(s2);
    q = -(a*c + 2).*c;
    dq = (u - acos(c))./s + (c - t*s2);
    q = q + dq;
    q = max(min(q,hi),lo);
end % Iteration

c = (sqrt(1- a*q)-1)/a;
s = sqrt(1-c.^2);
rho = R * (1 + a*c);
theta = (2*pi)*v(2,:);
% theta = pi/4 + (pi/4)*v(2,:);
x = cos(theta).*rho;
y = sin(theta).*rho;
z = r*s.*(2*v(3,:)-1);

% xyz = [x; y; z];s

end % torusmap 