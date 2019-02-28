# Matlab Reference

Quick reference for Matlab programming language.

Table of Contents
=================

   * [Matlab Reference](#matlab-reference)
   * [Table of Contents](#table-of-contents)
   * [Getting Help](#getting-help)
   * [Arrays](#arrays)
      * [Creating](#creating)
      * [Slicing](#slicing)
      * [Operations](#operations)
   * [User Interaction](#user-interaction)
      * [Output](#output)
      * [Input](#input)
   * [Loops](#loops)
      * [For Loop](#for-loop)
      * [While Loop](#while-loop)
   * [Conditional Statements](#conditional-statements)
      * [If Statement](#if-statement)
      * [Switch Statement](#switch-statement)
      * [Logic operators](#logic-operators)
   * [Functions](#functions)
   * [Plotting](#plotting)
      * [Single Figure](#single-figure)
      * [Multiple Figures](#multiple-figures)
      * [Line Types](#line-types)
      * [Subplot](#subplot)
   * [Control Systems](#control-systems)
      * [Transfer Functions](#transfer-functions)
         * [Defining](#defining)
         * [Analysing](#analysing)
      * [Block Diagrams](#block-diagrams)
   * [Reference](#reference)
      * [Mathematical Operations](#mathematical-operations)
      * [Data Analysis](#data-analysis)
      * [Time and Date](#time-and-date)
      * [Two-Dimensional Plotting](#two-dimensional-plotting)
      * [Three-Dimensional Plotting](#three-dimensional-plotting)
      * [Plot Annotation](#plot-annotation)
      * [Surface, Mesh and Contour Plots](#surface-mesh-and-contour-plots)
      * [Domain Generation and Interpolation](#domain-generation-and-interpolation)
      * [Polynomial Functions](#polynomial-functions)
      * [Discrete Time Signal Processing](#discrete-time-signal-processing)


# Getting Help

```Matlab
help command                        % Where command is what to find help for.
```
# Arrays

Arrays are indexed starting at 1, and referenced using parenthesis '()'.

## Creating

```Matlab
x = [1 4 6 8 3 5]                   % A list of numbers
x = 1:10                            % A range of numbers, start and end number included.
x = 1:0.1:5                         % A range of numbers with a set increment.
x = linspace(a,b,c)                 % A set of c numbers, between b and c (inclusive).
A = zeros(3,2)                      % 3x2 matrix of ones.
B = ones(4)                         % 4x4 matrix of zeros
R = rand(3)                         % 3x3 matrix of random numbers between 0 and 1
```

## Slicing

```Matlab
a = R(:,3)                          % All row column 3
b = R(2,:)                          % Row 2 all columns
c = R(2:3,[1 3])                    % Rows 2 and 3, columns 1 and 3
d = x(2:end)                        % Single dimension from index 2 to end 
```
## Operations

```Matlab
y = x.^2                            % Element wise operation require a period.
C = inv(R)                          % Inverse of a square matrix.
D = det(R)                          % Determinant of a square matrix. 
E = D`                              % Transpose matrix D.
F = isequal(E,D)                    % Check if two matrices are equal (0 or 1)
G = length(D)                       % Returns the length of the largest dimension.
H = size(E)                         % Return a vector of matrix dimension sizes (Rows,Columns..)

```        

# User Interaction

## Output

```Matlab
disp('Hello World')                 % Print a simple message.
disp(['joined' 'strings'])          % Strings can be joined if square brackes are used.
fprintf('A float. %3.2f\n',4.335)   % Formatted output.
int2str(123)                        % Convert integer to string.
```

## Input

```Matlab
in = input('Enter a number')        % Get input from user.
```

# Loops

## For Loop

```Matlab
for num = 1:10

    disp(int2str(num));
    for num = 1:10
        disp(int2str(num*2));
    end
end
```

## While Loop

```Matlab
vari = 0
while vari < 10

    disp(int2str(vari));
    vari = vari + 1;
end
```

# Conditional Statements

## If Statement

```Matlab
num = 5
if num < 2
    disp('Num less than 7')
elseif num < 7
    disp('Num less than 7')
else
    disp('Num greater than all')
end
```

## Switch Statement

```Matlab
num = 3
switch num
    case 1
        disp('Num = 1')
    case 2
        disp('Num = 2')
    case 3
        disp('Num = 3')
    otherwise
        disp('Num is huge')
end
```

## Logic operators

```
1 | 0                               % OR operator   (True)
1 & 0                               % AND operator  (False)
1 < 0                               % Less than     (False)
0 < 1                               % Greater than  (False)
0 == 1                              % Equal to      (False)
```

#  Functions

An example function, located in a separate file.

```Matlab
% FILE: myfunc.m
% File name must be the same as function name.
% File needs to be located in the same directory as script, or located on path.

function [r1 r2] = myfunc(a,b,c)    % Function takes 3 arguments and return 2 values.

r1 = a*b*c;                         % Colons used to supress output.
r2 = a+b+c;                         % Colons used to supress output.
```
Call the function from the main script like so:
```
[r1 r2] = myfunc(1,2,3)             % Perform the above operations on 1,2,3
```
# Plotting

## Single Figure

```Matlab
x = 1:10                            % Create an array
y = x.^2                            % Square the array
z = x.^3                            % Cube the array

figure                              % Make a new figure window
hold on                             % Freeze ploting so multiple line can be made.
plot(x,y)                           % Plot the first function
plot(x,z)                           % Plot the next function.
xlabel('x')                         % X-axis label.
xlabel('y')                         % Y-axis label.
title('A Title')                    % Plot title.
hold off                            % Plot everything.
```

## Multiple Figures

```Matlab
x = 1:10                            % Create an array.
y = x.^2                            % Square the array.
z = x.^3                            % Cube the array.

figure                              % Make a new figure window.
plot(x,y)                           % Plot the first figure.

figure
plot(x,z)                           % Plot the second figure.
```

## Line Types

Full list of line types on [Mathworks website](https://www.mathworks.com/help/matlab/ref/linespec.html).

```Matlab
x = 1:10                            % Create an array.
y = x.^2                            % Square the array.
z = x.^3                            % Cube the array.

figure                              % Make a new figure window.
hold on                             % Freeze ploting so multiple line can be made.
plot(x,y,'-')                       % Plot with a '-' line
plot(x,z,'+')                       % Plot with a '+' symbol
plot(x,z,':')                       % Plot with a dotted line.
hold off
```
## Subplot

```Matlab
x = 1:10                            % Create an array.
y = x.^2                            % Square the array.
z = x.^3                            % Cube the array.

figure                              % Make a new figure window
subplot(1,2,1)                      % 1x2 Subplot, position 1
plot(x,y)                           % Plot the first figure.
title('Plot 1')                     % Plot 1 title

subplot(1,2,2)                      % 1x2 Subplot, position 2
plot(x,z)                           % Plot the second figure.
title('Plot 2')                     % Plot 2 title
```

# Control Systems

## Transfer Functions

### Defining

```
% Method 1

numF1 = 24*[1,2,3];
denF1 = [1,0,4,0];
F1 = tf(numF1,denF1)

% Method 2

numF2 = poly([0 -2]);
denF2 = poly([-4 2 -5]);
F2 = tf(numF2, denF2)

% Method 3

s = tf('s');
F3 = s*(s + 2) / ((s+4) * (s-2) * (s+5))

```

### Analysing

```Matlab

Ts = s*(s + 2) / ((s+4)*(s-2))      % An arbritrary transfer function.

pole(Ts)                            % Find the poles of the transfer function.
zero(Ts)                            % Find the zeros of the transfer function.
pzmap(Ts)                           % Plot the poles and zeros in cartesian coordinates.
```

## Block Diagrams

```Matlab

G1s = tf([1],[2,3])                 % Forward path transfer function one.
G2s = tf([1],[2,3])                 % Forward path transfer function two.
H1 = tf([1],[1,0])                  % Feedback path tranfer function.

Gfp = series(G1s,G2s)               % Series combination, equivalent to G1s * G2s
Gfs = parallel(G1s,G2s)             % Parallel combination, equivalent to G1s + G2s
Ts = feedback(G1s,H1,-1)            % Combine feedback loop, integer represents positive/negitive.
```

# Reference

## Mathematical Operations

```Matlab
abs                                 % Absolute value complex magnitude.
acos, acosh                         % Inverse cosine, inverse hyperbolic cosine.
acot, acoth                         % Inverse cotangent, inverse hyperbolic cotangent.
acsc, acsch                         % Inverse cosecant,inverse hyperbolic cosecant.
angle                               % Phase angle.
asec, asech                         % Inverse secant, inverse hyperbolic secant.
asin, asinh                         % Inverse sine, inverse hyperbolic sine.
atan, atanh                         % Inverse tangent, inverse hyperbolic tangent
atan2                               % Four-quadrant inverse tangent.
ceil                                % Round toward infinity.
complex                             % Construct complex data from real and imaginary components.
conj                                % Complex conjugate.
cos, cosh                           % Cosine, hyperbolic cosine.
cot, coth                           % Cotangent, hyperbolic cotangent.
exp                                 % The exponential.
fix                                 % Round towards zero.
floor                               % Round towards minus infinity.
gcd                                 % Greatest common divisor.
imag                                % Imaginary part of a complex number.
lcm                                 % Least common multiple.
log                                 % Natural logarithm.
log2                                % Base 2 logarithm, dissect floating-point numbers.
log10                               % Common (base 10) logarithm.
mod                                 % Modulus.
nchoosek                            % Binomial coefficient or all combinations.
real                                % Real part of complex number.
rem                                 % Remainder after division.
round                               % Round to nearest integer.
sec, sech                           % Secant, and hyperbolic secant.
sign                                % Signum function, get sign of number.
sin, sinh                           % Sine, hyperbolic sine.
sqrt                                % Square root.
tan, tanh                           % Tangent, hyperbolic tangent.
```

## Data Analysis

```Matlab
corrcoef                            % Correlation coefficients.
cov                                 % Covarient matrix. 
cross                               % Vecotr cross product
cumprod                             % Cumulative product
cumsum                              % Cumulative sum
cumtrapz                            % Cumulative trapezoidal numerical integration
delaunay                            % Delaunay triangulation
del2                                % Discrete Laplacian
diff                                % Differences and approximate derivatives
fdsearch                            % Search for nearest point
factor                              % Prime factors
gradient                            % Numerical gradient
histogram                           % Histogram or Bar chart
inpolygon                           % Detect points inside a polygonal region
max                                 % Maximum elements of an array
mean                                % Average or mean value of arrays
median                              % Median value of arrays
min                                 % Minimum elements of an array
perms                               % All possible permutations
polyarea                            % Area of polygon
primes                              % Generate list of prime numbers
prod                                % Product of array elements
rand                                % Uniformly distributed random number
randn                               % Normally distributed random number
sort                                % Sort elements in ascending order
sortrows                            % Sort rows in ascending order
std                                 % Standard deviation
sum                                 % Sum of array elements
trapz                               % Trapezoidal numerical integration
tsearch                             % Search for enclosing Delaunay triangle
var                                 % Variance
voronoi                             % Voronoi diagram
```

## Time and Date

```Matlab
calendar                            % Calendar.
clock                               % Current time as a date vector.
cputime                             % Elapsed CPU time.
date                                % Current date string.
datenum                             % Serial date number.
datestr                             % Date string format.
datevec                             % Date components.
eomday                              % End of month.
etime                               % Elapsed time.
now                                 % Current date and time.
tic, toc                            % Stopwatch timer.
weekday                             % Day of the week.
```

## Two-Dimensional Plotting

```Matlab
bar                                 % Vertical bar chart.
barh                                % Horizontal bar chart.
hist                                % Plot histograms.
hold                                % Hold current graph.
loglog                              % Plot using log-log scales.
pie                                 % Pie plot.
plot                                % Plot vectors or matrices..
polar                               % Polar coordinate plot.
semilogx                            % Semi-log scale plot (X).
semilogy                            % Semi-log sclale plot (Y).
subplot                             % Crease multiple plots in one.
```

## Three-Dimensional Plotting        

```Matlab
bar3h                               % Vertical 3-D bar chart.
comet3                              % Horizontal 3-D bar chart.
cylinder                            % 3-D comet plot.
fill3                               % Generate cylinder.
plot3                               % Draw filled 3-D polygons in 3-space.
quiver3                             % Plot lines and points in 3-D space.
slice                               % 3-D quiver (or velocity) plot.
sphere                              % Volumetric slice plot.
stem3                               % Generate sphere.
waterfall                           % Plot discrete surface data.
```

## Plot Annotation

```Matlab
clabel                              % Add contour labels to a contour plot.
datetick                            % Date formatted tick labels.
grid                                % Grid lines for 2-D and 3-D plots.
gtext                               % Place text on a 2-D graph using a mouse.
legend                              % Graph legend for lines and patches.
plotyy                              % Plot graphs with Y tick labels on the left and right.
title                               % Titles for 2-D and 3-D plots.
xlabel                              % X-axis labels for 2-D and 3-D plots.
ylabel                              % Y-axis labels for 2-D and 3-D plots.
zlabel                              % Z-axis labels for 3-D plots.
```

## Surface, Mesh and Contour Plots

```Matlab
contour                             % Contour (level curves) plot.
contourc                            % Contour computation.
contourf                            % Filled contour plot.
hidden                              % Mesh hidden line removal mode.
meshc                               % Combination mesh/contourplot.
mesh                                % 3-D mesh with reference plane.
peaks                               % A sample function of two variables.
surf                                % 3-D shaded surface graph.
surface                             % Create surface low-level objects.
surfc                               % Combination surf/contourplot.
surfl                               % 3-D shaded surface with lighting.
trimesh                             % Triangular mesh plot.
trisurf                             % Triangular surface plot.
```

## Domain Generation and Interpolation

```Matlab
griddata                            % Data gridding and surface fitting.
meshgrid                            % Generation of X and Y arrays for 3-D plots.
polyfit(x,y,n)                      % Finds the coefficients of a polynomial of degreen n with x and y. 
interp1(x,Y,xi)                     % Returns vector elements corresponding to elements xi.
spline(x,y,xx)                      % Returns the f(xx) value at xx for interpolating cubic spline.
```

## Polynomial Functions

```Matlab
polyval(p1,x)                       % Polynomical (p1), evaluated at x.
conv(u,v)                           % Convolution.
deconv(v,u)                         % Deconvolution.
roots(c)                            % Returns a columns vector for roots of polynominal.
```

## Discrete Time Signal Processing

```Matlab
fft(xt)                             % Get the fast Fourier transform of x(t).
fftshift(xf)                        % Center a FFT about the zero frequency.
cheby2(A,B,C,D)                     % Get filter coeffecients for chebyshev filter.
filter(B,A,X)                       % Filter a time domain signal with filter coefficients.

```
