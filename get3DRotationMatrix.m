function R = get3DRotationMatrix(theta,dimension)
% returns the 3D rotation matrix

switch dimension
    case 1 % x rotation
        R = [1 0           0          ;...
             0 cos(theta) -sin(theta);...
             0 sin(theta)  cos(theta)];
    case 2 % y rotation
        R = [cos(theta) 0   sin(theta);...
             0          1   0         ;...
            -sin(theta) 0   cos(theta)];
    case 3 % z rotation
        R = [cos(theta) -sin(theta) 0;...
             sin(theta)  cos(theta) 0;...
             0           0          1];          

end