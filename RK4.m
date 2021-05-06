function [MY, MF] = RK4

load('params.mat')

MY = [];
MY(:, 1) = Y0;
MF = [];

for i = 1:length(t)-1
    
    Y = MY(:,i);
    
    K1 = [
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
        ];
    
    Y = MY(:,i) + h/2 * K1;
    
    K2 = [
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
        ];
    
    Y = MY(:,i) + h/2 * K2;
    
    K3 = [
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
        ];
    
    Y = MY(:,i) + h * K3;
    
    K4 = [
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
        ];
    
    Y = MY(:,i) + h/6 * (K1 + 2*K2 + 2*K3 + K4);
    MY(:, i+ 1) = Y;
    MF(:, i) = 1/6 * (K1 + 2*K2 + 2*K3 + K4);
    
end
end