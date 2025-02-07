%basic matrix operation 
inv(A)   inv(B)
A*B
A*A-B*B
A*A+2*A*B+B*B
Inv(A*B)
det(A)
det (B-A)* inv(A-B)*

%plot the vector field
x=linspace(-3,3,10);
y=linspace(-3,3,10);
z=linspace(-3,3,10);
[X,Y,Z]=meshgrid(x,y,z);
f1=11*Y.*Z;
f2= (3*X.^2) .*Z;
f3=5*X.^3;
quiver3(X,Y,Z,f1,f2,f3)

%fit a linear equation for the data
x=[0 1 2 3 4]
y=[5 7 8 11 12]
p=polyfit(x,y,1);
y_fit=polyval(p,x);
plot(x,y,'ro')
hold on
plot(x,y_fit,'b')

%Plot the standard functions:
%Plot x^2  in (0,10)
x=linspace(0,10,100);
y=x.^2 ;
y1=exp(x);
plot(x,y,'r')
hold on
plot(x,y1,'g')

%Plot the graph of y1=sinx and y2= cosx in (0,4)
x=linspace(0,4);
y=sin(x);
y1=cos(x);
plot(x,y,x,y1)
xlabel('x-axis');
ylabel('y-axis');
title('Graph of sin(x) and cos(x)');
legend('sin(x)','cos(x)');

%Numerical  Integration - Trapezoidal 
function [  ] = Trep_int(x0,xn,n)
h=(xn-x0)/n;
f=@(x) 1./(1+x.^2);
sum=f(x0)+f(xn);
for i=1:n-1
    x=x0+i*h;
    sum=sum+2*f(x);
end
I=(h/2)*sum;
disp(I)
end

%Numerical Integration - Simpson 1/3rd
function [  ] = Simpsons1_3rule(x0,xn,n)
h=(xn-x0)/n;
f=@(x) 1./(1+x.^2)
sum=f(x0)+f(xn);
for i=1:n-1
    x=x0+i*h;
    if mod(i,2)==0
        sum=sum+2*f(x);
    else
        sum=sum+4*f(x);
    end
end
I=(h/3)*sum;
disp(I)
end

%Numerical Interpretation - Simpson 3/8th
function [  ] = Simpsons3_8rule(x0,xn,n)
h=(xn-x0)/n;
f=@(x) 1./(1+x.^2); 
sum=f(x0)+f(xn);
for i=1:n-1
    x=x0+i*h;
    if mod(i,3)==0
        sum=sum+2*f(x);
    else
        sum=sum+3*f(x);
    end
end
I=(3*h/8)*sum;
disp(I)
end

%gauss-seidel
A = [2 1 1; 3 5 2; 2 1 4];
B = [5; 15; 8]; 
 
X = zeros(size(B)); 
n = length(B); 
tolerance = 1e-6; 
max_iterations = 100;
 
for k = 1:max_iterations 
    X_old = X;  
    for i = 1:n 
        sum1 = A(i, 1:i-1) * X(1:i-1);
        sum2 = A(i, i+1:n) * X_old(i+1:n);
         
        X(i) = (B(i) - sum1 - sum2) / A(i, i); 
    end 
     
    if norm(X - X_old, inf) < tolerance 
        fprintf('Converged in %d iterations.\n', k); 
        break; 
    end 
end 
 
if k == max_iterations 
    fprintf('Did not converge within the maximum number of iterations.\n'); 
else 
    fprintf('Solution:\n'); 
    disp(X); 
end 

%lagranze interpolation 
function [poly] = lag_poly(x,y) 
n=length(x); 
syms t 
poly=0; 
for i=1:n 
    L=1; 
    for j=1:n 
        if(i~=j) 
        end 
    end 
    L=L*(t-x(j))/(x(i)-x(j)); 
    poly=poly+L*y(i); 
end 
poly=simplify(poly); 
t = x(1):0.001:max(x); 
z = eval(poly); 
%plot(x,y,'*',t,z) 
plot(x,y,'*k',t,z,'k') 
legend('Data Points','Lagrange Polynomial') 
xlabel('x') 
ylabel('y,p(x)') 
grid on 
end 

%Newton Raphson
function NewtonRaphson(f, df, x0)
x(1) = x0;
tolerance = 10^-6;
for i = 1:100
    if abs(df(x(i))) < 10^-15
        fprintf('Derivative of f is zero')
        break;
    end
else
    x(i+1) = x(i) - f(x(i))/df(x(i));
    if abs(x(i+1) - x(i)) < tolerance
        fprintf('The root of the given equation is %f \n', x(i));
        break;
    end
end
end


