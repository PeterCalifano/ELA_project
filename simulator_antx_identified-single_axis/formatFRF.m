function FRF_ReIm = formatFRF(FRF)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Determine size of input
[nrow, ncol] = size(FRF);
% Transpose if frequency axis is along the rows
if nrow < ncol
    FRF = FRF';
    [nrow, ncol] = size(FRF);
end

FRF_ReIm = zeros(2*nrow, ncol);
% Loop over the entries and split Re and Im part
for i = 1:ncol
    FRF_ReIm(:, i) = [real(FRF(:, i)); imag(FRF(:, i))];
end

end