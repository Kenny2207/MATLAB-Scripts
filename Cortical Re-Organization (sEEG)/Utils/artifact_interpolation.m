function data_3D = interpolate_artifact_hardcoded(data_3D, artifact_indices, valid_window_size)
% INTERPOLATE_ARTIFACT_HARDCODED - Performs localized interpolation of stimulation artifacts.
%
% INPUTS:
%   data_3D          - 3D Matrix containing neural data (channels x time-samples x trials).
%   artifact_indices - Vector of time-bin indices targeted for interpolation (e.g., 343:347).
%   valid_window_size- (OPTIONAL) Number of clean samples to consider before and after 
%                      the artifact window for boundary estimation. Default: 10.
%
% OUTPUT:
%   data_3D          - 3D Matrix containing the artifact-corrected timeseries data.
    
    % --- Assign default parameters if missing ---
    if nargin < 3 || isempty(valid_window_size)
        valid_window_size = 10;
    end
    
    % --- Determine boundaries for localized interpolation ---
    indices_before = (artifact_indices(1) - valid_window_size) : (artifact_indices(1) - 1);
    indices_after  = (artifact_indices(end) + 1) : (artifact_indices(end) + valid_window_size);
    valid_indices  = [indices_before, indices_after];
    
    [n_channels, ~, n_trials] = size(data_3D);
    
    % --- Execution Loop across Trials and Channels ---
    for i_trial = 1:n_trials
        for i_chan = 1:n_channels
            % Extract the raw timeseries trace for the current channel/trial configuration
            signal_trace = squeeze(data_3D(i_chan, :, i_trial));
            
            % Execute Piecewise Cubic Hermite Interpolating Polynomial ('pchip') 
            % to ensure shape-preserving, non-oscillatory waveform reconstruction
            interpolated_values = interp1(valid_indices, signal_trace(valid_indices), artifact_indices, 'pchip');
            
            % Substitute the artifactual data segment with the interpolated waveform
            data_3D(i_chan, artifact_indices, i_trial) = interpolated_values;
        end
    end
    
    fprintf('Artifact correction successfully executed on indices %d:%d.\n', artifact_indices(1), artifact_indices(end));
end