function [MY, MF] = Euler()

load('params.mat')

MY = [];
MY(:, 1) = Y0;
MF = [];

for i=1:length(t)-1
    
    Y = MY(:,i);
    
    K = [
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
        ];
    
    Y = MY(:,i) + h * K;
    
    MY(:,i+1) = Y;
    MF(:,i) = K;
    
end
end