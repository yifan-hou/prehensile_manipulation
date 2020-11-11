function [R] = rotz_deg(deg)
ang = deg*pi/180;
c = cos(ang);
s = sin(ang);
R = [c -s 0;
     s c  0;
     0 0  1];
end

