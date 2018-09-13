% Project 2
% Karol Jose Gutierrez Suarez A01024536
% Raziel Nicolas Martinez Castillo A01410695
% Emilio Hernández López A01336418
%Jesus Bernardo Granados 
load('SatuSignal.mat');
x=1;

%declare a matrix for data in X and Y such that Y=f(X)
X=tsat;
Y=fsat;

%the data consists of 275 points
size(X);%275
size(Y);%275

%This matrix stores the roots
roots=[];

%roots are obtained using biection method taking the closest two points to the
%X axis such that they are in different sides of the X axis
for i=2:275
       if Y(i).*Y(i-1)<=0%this way we know that the points are in different sides of the X axis
           roots=[roots, X(i)-(Y(i).*(X(i)-X(i-1))./(Y(i)-Y(i-1)))];
       end
end

%plot(roots, 0, "^"), hold on
roots;
valuec=[];
size(roots); %21
%We can choose two consecutive roots to solve the c value, this comes from
%the fact that the sine function has a value of zero every PI units.
for i=2:21
    valuec=[valuec, pi./((roots(i).^2-roots(i-1).^2))];
end
%the c values are a different but similar due to the bisection, so we take
%the average
global C;
C=mean(valuec);
%The next step is to solve A and B. 

%We consider the crests and valleys because these values are local
%extreme values and then it must happen that sin(cx^2)=+-1 and then the
%function in these points become f(x)=+-Ae^(-Bx)
%Crests and valleys are stored in the matrix points
points=[];

%Fill the matrix
for i=2:274
   if Y(i)>Y(i-1)&&Y(i)>Y(i+1)%case crest
       points=[points, i];
   end
   if Y(i)<Y(i-1)&&Y(i)<Y(i+1)%case trough
       points=[points, i];
   end
end

%plot(x(points),y(points), "^"), hold on

%use exponential model regresion to find the curve Ae^(-Bx) using the
%minimum and maximum points

%points
%as in the formula for the exponential model, we have vaiables z and w
z=[];
w=X(points);%according to the formula, w value is equal to the x value of the function
size(points);%11

for i=1:11%obtain z value using log
    z= [z, log(  Y(points(i))  ./  sin(C* X(points(i)).^2)) ];%we divide over the sin() value to adjust the size and orientation of the function
end
%z
%w
%here we store the result of the matrix multiplication that lead to the
%values b0 and b1
b0b1=[sum(z) sum(w.*z)]*[11 sum(w); sum(w) sum(w.^2)]^(-1);

%from the formula for exponential model, we know which are the values for A
%and B
global A;
A=exp(b0b1(1));
global B;
B=-b0b1(2);%negative because our function has Ae^(-Bx)

A;
B;
C;
f=A.*exp(-B.*x).*sin(C.* x.^2);

%plot the points
%plot(tsat, fsat, "-"), hold on

figure('Name','Original function')
plot(tsat, fsat, "*"), hold on
dx=0.001;
x=0: dx:3.5;
plot(x, f), hold on

fd=A.*exp(-B .*x).*(2 *C.*x.*cos(C .*x.^2) - B.*sin(C.*x.^2));
figure('Name','First derivative')
plot(x, fd), hold on

f2d=A.*exp(-B.*x).*((B.^2 - 4*C.^2.*x.^2).*sin(C.*x.^2) - 2.*C.*(2.*B.*x - 1).*cos(C.*x.^2));
figure('Name','Second derivative')
plot(x, f2d), hold on
%we can obtain the values of the crests and valleys 
fprintf('\n\nCrests and valleys are:\n');
for i=1:7
    if mod(i, 2)==1
        fprintf('Crest:\n  x=%f\n   y=%f\n\n', sqrt((i-0.5)*pi/C), f1(sqrt((i-0.5)*pi/C)));
    else
         fprintf('Valley:\n  x=%f\n   y=%f\n\n', sqrt((i-0.5)*pi/C), f1(sqrt((i-0.5)*pi/C)));
    end
end

%inflexion points
infPoints=[];
i=0;
%we use the second derivative and find the roots of this function
while length(infPoints)<=4 && i<=3 
    %we use bisection method to approximate the root values
    while(f3(i)*f3(i+0.05)>0)
        i=i+0.05;
    end
    xi=i+0.05; 
    xim=i;
    xn=(xi+xim)/2;
    for j=1:50
        if f3(xn).*f3(xi) < 0
            xim=xn;
        elseif f3(xn).*f3(xi) < 0
                xi=xn;
        end
        xn=(xi+xim)./2;    
    end
    %eliminate duplicates
    infPoints=notRepeated(infPoints, xn, 0.0001);
    i=i+0.05;
end
%print results
fprintf('\n\nFirst inflexion points are:\n');
for i=1:5
    fprintf('Inflexion point %d  x=%f  y=%f\n\n', i, infPoints(1, i), f1(infPoints(1, i))); 
end
%function for ptoblem 4
fp4=A.*exp(-B.*x).*sin(C.* x.^2)-1.5;
dx=0.001;
x=0: dx:3.5;
figure('Name','Function for P4')
plot(x, fp4), hold on
arrp4=[];
i=0;
%we procede using bisection method to find the roots of the new function
while i<=3 
    while f4(i)*f4(i+0.05)>0 && i<=3 
        i=i+0.05;
    end
    xi=i+0.05; 
    xim=i;
    xn=(xi+xim)/2;
    for j=1:50
        if f4(xn).*f4(xi) < 0
            xim=xn;
        elseif f4(xn).*f4(xi) < 0
                xi=xn;
        end
        xn=(xi+xim)./2;    
    end
   
    arrp4=notRepeated(arrp4, xn, 0.0001);
    i=i+0.05;
end
arrp4;
sump4=0;
%using the graph we know that we need to sum the length ogf all the
%intervals such that the end of the interval has an even number as index
i=2;
while i<length(arrp4)
    sump4=sump4+arrp4(i)-arrp4(i-1);
    i=i+2;
end
sump4;
%print results
fprintf('\n\nThe when the signal was above 1.5: %f\n', sump4);

%fucntion for problem 5
fp5=A.*exp(-B.*x).*sin(C.* x.^2)+1;
dx=0.001;
x=0: dx:3.5;
figure('Name','Function for P5')
plot(x, fp5), hold on

arrp5=[];
i=0;
%we procede using bisection method to find the roots of the new function
while i<=3 
    
    while f5(i)*f5(i+0.05)>0 && i<=3 
        i=i+0.05;
    end
    xi=i+0.05; 
    xim=i;
    xn=(xi+xim)/2;
    for j=1:50
        if f5(xn).*f5(xi) < 0
            xim=xn;
        elseif f5(xn).*f5(xi) < 0
                xi=xn;
        end
        xn=(xi+xim)./2;    
    end
   
    arrp5=notRepeated(arrp5, xn, 0.0001);
    i=i+0.05;
end
arrp5;
sump5=0;
i=2;
%using the graph we know that we need to sum the length ogf all the
%intervals such that the end of the interval has an even number as index
while i<length(arrp5)
    sump5=sump5+arrp5(i)-arrp5(i-1);
    i=i+2;
end
sump5;
%print results
fprintf('\n\nThe total time when the signal was below -1: %f\n', sump5);

%problem 6
i=1;
extremes=[];
%we store the values of crests and troughs using the fact that
%sin(cx^2)=+-1 in these points so we can solve for x
while i<50
    extremes=[extremes, sqrt((i-0.5)*pi./C)];
    i=i+1;
end
i=1;
%we find the first value that satisfies the desired condition
while abs(f1(extremes(i)))>0.1
i=i+1;
end

fprintf('\n\nThe time when the signal amplitude is below 0.1 is: %f\n', extremes(i));


function arr=notRepeated(arr, x, tol)
    repeated=false;
    for i=1:size(arr, 2)
        if abs( arr(i)-x)<tol
            repeated=true;
            break
        end    
    end
    if ~repeated
       arr=[arr, x]; 
    end
end
 
%function f
function y=f1(x)
    global A B C;
    y=A.*exp(-B.*x).*sin(C.* x.^2);
end


%second derivative
function y=f3(x)
    global A B C;
    y=A.*exp(-B.*x).*((B.^2 - 4*C.^2.*x.^2).*sin(C.*x.^2) - 2.*C.*(2.*B.*x - 1).*cos(C.*x.^2));
end

%function f4
function y=f4(x)
    global A B C;
    y=A.*exp(-B.*x).*sin(C.* x.^2)-1.5;
end

%function f5
function y=f5(x)
    global A B C;
    y=A.*exp(-B.*x).*sin(C.* x.^2)+1;
end


