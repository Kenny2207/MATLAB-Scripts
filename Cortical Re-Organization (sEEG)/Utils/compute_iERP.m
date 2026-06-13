function [erp_data, res_leads_ierp, v_times_ierp] = compute_iERP(data, cnames, srate, low_cutoff, high_cutoff, order, pthr, baseline_range, post_baseline_range, apply_interpolation, idx_artifact, idx_valid)
% COMPUTE_IERP: Computes intracranial Evoked Responses (iERP) and performs cluster-based significance testing.
%
% PIPELINE WORKFLOW:
%   1. Optional localized interpolation for artifact removal (e.g., stimulation artifacts).
%   2. Zero-phase bandpass filtering (Butterworth filter).
%   3. Downsampling to reduce computational overhead and align timescales.
%   4. Mean-subtraction baseline correction per individual trial.
%   5. Trial-by-trial Z-score normalization relative to the baseline period variance.
%   6. Grand average ERP computation across trials and one-sample t-test (vs. 0).
%   7. Temporal clustering to identify channels with sustained significant activations.
%
% INPUT PARAMETERS:
%   data                - 3D numeric matrix (channels x time-samples x trials) containing raw sEEG/EEG data.
%   cnames              - Cell array of strings containing channel labels.
%   srate               - Integer, current sampling rate in Hz (e.g., 1000 Hz).
%   low_cutoff          - High-pass filter cutoff frequency (Hz).
%   high_cutoff         - Low-pass filter cutoff frequency (Hz).
%   order               - Order of the Butterworth filter.
%   pthr                - Alpha significance threshold for statistical testing (e.g., 0.05).
%   baseline_range      - Vector of indices defining the baseline epoch window (e.g., 1:50).
%   post_baseline_range - Vector of indices defining the post-stimulus window targeted for evaluation.
%   apply_interpolation - Boolean flag (true/false) to enable/disable artifact interpolation.
%   idx_artifact        - (Optional) Time indices targeted for interpolation. Required if apply_interpolation is true.
%   idx_valid           - (Optional) Clean boundary time indices used for interpolation anchors. Required if apply_interpolation is true.

    %% --- Filter Parameter Initialization ---
    [b, a] = butter(order, [low_cutoff, high_cutoff] / (srate / 2), 'bandpass');
    [n_channels, ~, n_trials] = size(data);
    
    % Constant added to baseline variance to prevent division-by-zero runtime errors
    EPSILON = 1e-6;

    % Output structure initialization
    erp_data = struct;
    res_leads_ierp = {};

    %% --- Timeline & Window Segment Definition ---
    downsample_factor = 5;
    v_times_ierp = linspace(-300, 699, 1000 / downsample_factor);

    % Segment boundaries to extract the primary 1000-sample core epoch (-300 ms to +699 ms)
    start_idx = 1; 
    end_idx = 1000;

    %% --- Channel Processing Loop ---
    for i = 1:n_channels
        % Extract the core 1000-sample segment [1000 samples x n_trials]
        data_signal = squeeze(data(i, start_idx:end_idx, :));  

        % --- 1. OPTIONAL ARTIFACT INTERPOLATION (Executed prior to filtering to prevent ringing artifacts) ---
        if apply_interpolation
            if nargin < 12
                error('Error: In order to execute interpolation, both idx_artifact and idx_valid must be provided.');
            end
            
            % Restrict artifact and anchor indices within the 1000-sample sub-window boundaries
            valid_points = idx_valid(idx_valid >= 1 & idx_valid <= size(data_signal, 1));
            artifact_points = idx_artifact(idx_artifact >= 1 & idx_artifact <= size(data_signal, 1));
            
            if length(valid_points) >= 2 && ~isempty(artifact_points)
                for trial = 1:n_trials
                    single_trial_trace = data_signal(:, trial);
                    % Linear extrapolation/interpolation to reconstruct contaminated segments
                    interpolated_values = interp1(valid_points, single_trial_trace(valid_points), artifact_points, 'linear', 'extrap');
                    data_signal(artifact_points, trial) = interpolated_values;
                end
            end
        end
        
        % --- 2. SIGNAL FILTERING & PRE-PROCESSING ---
        
        % Apply zero-phase forward and reverse Butterworth filter
        data_filt = filtfilt(b, a, data_signal);

        % Decimate timeseries via downsampling -> [200 samples x n_trials]
        data_ds = downsample(data_filt, downsample_factor);

        % Store channel label
        erp_data(i).label = cnames{i};

        % Baseline Correction: Subtract the trial-specific baseline mean
        baseline_mean_per_trial = mean(data_ds(baseline_range, :), 1);
        data_corr = data_ds - baseline_mean_per_trial;

        % Trial-by-Trial Z-score Normalization based on baseline variance
        data_z = zeros(size(data_corr));
        for trial = 1:n_trials
            baseline_vals = data_corr(baseline_range, trial);
            mu = mean(baseline_vals, 'omitnan');
            sigma = std(baseline_vals, 0, 'omitnan');
            
            % Normalize trace using the epsilon factor for numerical stability
            data_z(:, trial) = (data_corr(:, trial) - mu) / (sigma + EPSILON);
        end

        % Compute grand average ERP across all trials
        erp_avg = mean(data_z, 2, 'omitnan');
        erp_data(i).compute = erp_avg;
        erp_data(i).z = data_z; % Store normalized trials for downstream diagnostics

        %% --- 3. STATISTICAL EVALUATION ---
        erp_data(i).p = ones(1, length(post_baseline_range));
        erp_data(i).t = zeros(1, length(post_baseline_range));
        
        for idx = 1:length(post_baseline_range)
            tp = post_baseline_range(idx);
            signal_tp = data_z(tp, :);
            % One-sample t-test against zero across trials
            [~, erp_data(i).p(idx), ~, stats] = ttest(signal_tp, 0);
            erp_data(i).t(idx) = stats.tstat;
        end
        
        % Track all significant time bins (capturing both positive and negative deflections)
        erp_data(i).idxs = find(erp_data(i).p < pthr); 
        erp_data(i).comparison = erp_data(i).idxs;
        
        % Cluster consecutive significant time-bins to control for false positives
        if ~isempty(erp_data(i).comparison)
            erp_data(i).grouped = mat2cell(erp_data(i).comparison, 1, ...
                diff([0, find(diff(erp_data(i).comparison) ~= 1), length(erp_data(i).comparison)]));
            
            for g = 1:length(erp_data(i).grouped)
                % Channel is marked as responsive if a cluster contains at least 3 consecutive significant bins
                if length(erp_data(i).grouped{g}) >= 3
                    if ~ismember(erp_data(i).label, res_leads_ierp)
                        res_leads_ierp{end+1} = erp_data(i).label;
                    end
                end
            end
        else
            erp_data(i).grouped = {}; 
        end
    end

    % --- 4. ARTIFACT REJECTION & EXTRA-CRANIAL FILTERING ---
    excluded_channels = ["EOG", "EKG", "DC01", "DC02", "DC03", "DC04", "DC09", "DEL SX", "DEL DX"];
    res_leads_ierp = res_leads_ierp(~ismember(string(res_leads_ierp), excluded_channels));
end