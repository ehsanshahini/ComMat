function R = RotatM(v)
% R = RotatM(v)
% Given a 4D-vector with first element equal to rotational angle and rest
% thress element equal to rotational axis, this program generate
% Rotational matrix.

if numel(v)<4
    disp('the argument must be a 4D-vector!!')
    return
end
v=v(1:4);

v(2:4) = v(2:4)/norm(v(2:4));

a = cos(v(1)/2);
b = sin(v(1)/2)*v(2);
c = sin(v(1)/2)*v(3);
d = sin(v(1)/2)*v(4);

R = [a^2+b^2-c^2-d^2 2*b*c-2*a*d 2*a*c+2*b*d
     2*a*d+2*b*c a^2-b^2+c^2-d^2 2*c*d-2*a*b
     2*b*d-2*a*c 2*a*b+2*c*d a^2-b^2-c^2+d^2];