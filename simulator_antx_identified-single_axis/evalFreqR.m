function [TF_FRF] = evalFreqR(TF, freqs, unit)
%% PROTOTYPE
% [TF_FRF, TFmr] = evalFreqR(TF, freqs, unit)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate Frequency Response Function of Transfer Function Matrix TF (of
% any size for MIMO systems) over a grid of frequencies specified by freqs 
% expressed in "unit" measurement unit.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% TF: [NxM] Transfer Function MATLAB class with N rows and M columns
% freqs: [Nfx1] Frequency grid for evaluation 
% unit: [char] Measurement unit input to freqresp() function. Typical 'Hz'
%              or 'rad/s'. Default: 'rad/s'.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% TF_FRF: [(2*Nfx(N*M)] Frequency Response Function. Each column contains 
%  [Re(FRF); Im(FRF)] of each entry of the TF (indexing along columns as 
%  done by MATLAB).
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01/02/2023    Pietro Califano    First version - Coded and tested
%
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Control System Toolbox
% formatFRF() function
% -------------------------------------------------------------------------------------------------------------

%% Function code

if nargin < 3
    unit = 'rad/s';
end

% Apply minreal to TF to use minimal realization
% TFmr = minreal(TF);
TFmr = TF;

% Determine number of frequencies
Nf = length(freqs);
% Determine size of the TF
[nrow, ncol] = size(TFmr);
Ntf = nrow*ncol;

% Static allocation of output
% NOTE: Each column contains [Re(FRF); Im(FRF)] of each entry of the TF (indexing
% along columns as done by MATLAB
TF_FRF = zeros(2*Nf, Ntf);

for i = 1:Ntf
    % Evaluate frequency response over freqs grid
    TFi = squeeze(freqresp(TFmr(i), freqs, unit));
    % Allocate real and imag part of FRF in output vector
    TF_FRF(:, i) = formatFRF(TFi);

end

end