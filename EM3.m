%% matrix 
A = [ 9 2 1 ; 2 4 6 ; 1 5 7 ];

disp("(A^3") - 3*A");
disp((A^3) - 3*A);

disp("determinant of A");
disp(det(A)); 

disp("Adjacent of A if exists");
if det(A) ~=0;
    disp(det(A) * inv(A));
else
    disp("It's a singular matrix");
end

disp("Inverse of A");
disp(inv(A));

disp("Transpose of A");
disp(A');

%% plot the given function 
f = @(x) sin(x);                 %% acc to question
pts = linspace(0, 4*pi, 1000);
y = f(pts);
plot(pts,y,'b');

xlabel("X axis");
ylabel("Y axis");
title("Sin(x) graph");
grid on; 

%%plot the vector field 
x = linspace(-3, 3, 10);  %% change x y z acc to question 
y = linspace(-3, 3, 10);
z = linspace(-3, 3,10);
[x, y, z] = meshgrid(x, y, z);

u = x.* cos(x);
v = y.* sin(y);
w = z.*(z-y);

quiver3(x, y, z, u, v, w);

xlabel("X axis");
ylabel("Y axis");
zlabel("Z axis");
title("Vector plot graph");

%% polyfit (2a)
x = [10, 12, 16, 20];       %%change acc to question
y = [29, 33, 41, 53];

p = polyfit(x,y,2);
yfit = polyval(p,x);
plot(x,y,'ro');

hold on;
plot(x,yfit,'b');

xlabel("X axis");
ylabel("Y axis");
title("Fitting the curve"); 

%%lagrange's interpolation (a)
function [poly] = lag_pol(x,y)
n = length(x);
t_vals = linspace(min(x), max(x),1000);
poly = zeros(1, length(t_vals));
for i = 1:n
    L = ones(1, length(t_vals));

    for j = 1:n
        if i ~= j
            L = L.*(t_vals - x(j)) / (x(i) - x(j));
        end
    end
    poly = poly + L *y(i);
end
plot(x,y,'*k', t_vals, poly, 'k');
legend('Data points', Lagrange Polynomial');
xlabel('X - axis');
ylabel('Y - axis, p(x)');
grid on;
end

%%3(b) 
x = [0, 1, 2, 4];       %% change acc to question
y = [-1, 2, 7, 23];

p = polyfit(x, y, 2);   %% fitting 2nd degree polynomial
yfit = polyval(p,x);    %% evaluate the polynomial with x points
plot(x, y, 'ro');

hold on;
plot(x, yfit,'b');

xlabel("X axis");
ylabel("Y axis");
title("Fitting the curve");

f3 = polyval(p,3);
disp("Estimated value of f(3) : ");
disp(f3) 

%%simpson 1/3 
function [] = simp1_3(x0, xn, n)
h = (xn-x0)/n;
f = @(x) 1 ./ (1+x); %%change acc to question
sum = f(x0)+f(xn);
for i=1:n-1
    x = x0 + i * h;
    if mod(i,2)==0
        sum=sum+2*f(x);
    else
        sum=sum+4*f(x);
    end
end
I = (h/3)*sum;
disp(I)
end

%%taking 6 sub interval 
functon [] = intrvalOnly(x0, xn, n)
  h = (xn = x0) / n;
  f=@(x) 1 ./(1+x); %%change acc to question 
  sum = f(x0) + f(xn);
  for i=1:n-1
      x = x0 + i * h;
      sum sum + 2 * f(x);
  end
  I=(h/2)* sum;
  disp(I)
end 
