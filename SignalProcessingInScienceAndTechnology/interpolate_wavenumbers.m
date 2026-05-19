function [f_vec,k_a, k_s ] = interpolate_wavenumbers(exct,f_A0, k_A0,f_S0,k_S0  )
% interpolate_wavenumbers - Interpolate acoustic wavenumbers onto simulation frequency grid
%
% Syntax:
%   [f_vec, k_a, k_s] = interpolate_wavenumbers(exct, f_A0, k_A0, f_S0, k_S0)
%
% Description:
%   This function prepares frequency and wavenumber vectors for simulation by
%   constructing the simulation frequency axis and interpolating given dispersion
%   data (compressional and shear wavenumbers) onto that axis.
%
%   Inputs:
%     exct  - structure with simulation parameters:
%               exct.fs     : sampling frequency in Hz
%               exct.n_sampl: number of frequency samples (positive integer)
%     f_A0  - vector of frequencies (kHz) corresponding to compressional wavenumbers
%     k_A0  - vector of compressional wavenumbers (1/m) corresponding to f_A0
%     f_S0  - vector of frequencies (kHz) corresponding to shear wavenumbers
%     k_S0  - vector of shear wavenumbers (1/m) corresponding to f_S0
%
%   Outputs:
%     f_vec - column vector of simulation frequencies in Hz (0 to fs-df)
%     k_a   - column vector of compressional wavenumbers interpolated onto f_vec (rad/m)
%     k_s   - column vector of shear wavenumbers interpolated onto f_vec (rad/m)
%
%   Notes:
%     - Input frequency vectors f_A0 and f_S0 are expected in kilohertz (kHz)
%       and are converted to Hz inside the function.
%     - Input wavenumbers k_A0 and k_S0 are expected in 1/m and are converted
%       to rad/m (multiplied by 1e3 when matching units in subsequent code).
%     - Interpolation uses shape-preserving piecewise cubic interpolation ('pchip').
%       Values outside the provided dispersion ranges are initially set to NaN and
%       then replaced by the nearest available value to avoid gaps.
%
%   Examples:
%     % See calling code for typical usage where exct.fs and exct.n_sampl are set.
%
%   Author: Generated helper description
%   Date:   2026
%
% Prepare output frequency vector (Hz) for the rest of the function
f_vec = (0:exct.n_sampl-1)' * (exct.fs*2 / exct.n_sampl);  % column vector 0..fs-df

% Preallocate outputs k_a and k_s to match f_vec size (filled later)
k_a = nan(size(f_vec));
k_s = nan(size(f_vec));
df = exct.fs*2 / exct.n_sampl;
f_sim = (0:exct.n_sampl-1)' * df;        % column vector, 0..fs-df

% Original dispersion data use f_vec (Hz) and k_a, k_s (rad/m) from loaded file.
% Interpolate k vectors onto f_sim. Use NaN for out-of-range values then fill with nearest.
k_a_sim = interp1(f_A0*1e3, k_A0*1e3, f_sim, 'pchip', NaN);
k_s_sim = interp1(f_S0*1e3, k_S0*1e3, f_sim, 'pchip', NaN);

% Fill NaNs at ends by nearest valid value to avoid gaps in further processing
if any(isnan(k_a_sim))
    valid = find(~isnan(k_a_sim));
    if ~isempty(valid)
        k_a_sim(1:valid(1)-1) = k_a_sim(valid(1));
        k_a_sim(valid(end)+1:end) = k_a_sim(valid(end));
    end
end
if any(isnan(k_s_sim))
    valid = find(~isnan(k_s_sim));
    if ~isempty(valid)
        k_s_sim(1:valid(1)-1) = k_s_sim(valid(1));
        k_s_sim(valid(end)+1:end) = k_s_sim(valid(end));
    end
end

% % Save interpolated vectors into variables expected later (f_vec, k_a, k_s)
% % Overwrite f_vec to the simulation frequency grid and ensure column orientation
f_vec = f_sim;
k_a = k_a_sim;
k_s = k_s_sim;
