function [X] = skew(x)
%   Performs the cross product or skew symetric form of a 3x1 vector.

X = [ 0    -x(3)  x(2);
      x(3)  0    -x(1);
     -x(2)  x(1)  0];

end