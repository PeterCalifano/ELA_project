function fcn = Hmodelstruct(th)
%% PROTOTYPE
% fcn = Hmodelstruct(th)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function returning the transfer functions of the Longitudinal dynamics
% given theta parameters vector
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% th: [6x1] parameters vector
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% fcn: [2x1] Transfer functions of the system
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     Function documented
% -------------------------------------------------------------------------------------------------------------


%% Function code
s = tf('s');
den = 9.81*th(3) - s^2 *th(1) - s^2 *th(4) + s^3 + s*th(1)*th(4) - s*th(2)*th(3);


fcn = [s*(th(3)*th(5) - th(6)*th(1) + th(6)*s)/den; 

       (9.81*th(3)*th(5) - 9.81*th(6)*th(1) + th(5)*s^3  + th(6)*th(2)*s  - th(4)*th(5)*s)/den];

end
