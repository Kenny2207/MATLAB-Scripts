function [tf_data, res_leads] = time_freq_analysis(data, cnames, s_event, srate, v_times, pthr, baseline_range, post_baseline_range, condition)
% TIME_FREQ_ANALYSIS: Executes time-frequency decomposition on an EEG/sEEG dataset.
%
% METHODOLOGICAL NOTE: Artifact interpolation should be performed on absolute power 
% (transformed data) PRIOR to Z-score normalization. This ensures that the baseline 
% variance estimation is not contaminated by high-amplitude artifactual transients.

    % --- Analysis Parameters ---
    MIN_CLUSTER_LENGTH = 3; % Minimum consecutive time bins to define a significant response
    EPSILON = 1e-6;         % Small constant to prevent division by zero during Z-scoring

    % --- Initialization ---
    n_channels = size(data,1);
    n_trials = size(data, 3);
    tf_data = struct;
    res_leads = {};

    for i=1:n_channels
        % Extract signal for the current channel
        data_power = squeeze(data(i,:,:));
        tf_data(i, 1).label = cnames{i}; % Cell indexing ({}) ensures robustness for string arrays
        
        % --- 1. SPECTRAL POWER DECOMPOSITION ---
        [tf_data(i,1).transform, tf_data(i,1).freqs_out, tf_data(i,1).times_out] = timefreq(data_power, srate, 'freqs', [55, 145], 'timesout', v_times, 'cycles', 4, 'nfreqs', 10);
        tf_data(i).transform_abs = abs(tf_data(i).transform);

        % --- 2. Z-SCORE NORMALIZATION (Performed on clean absolute power data) ---
        tf_data(i).norm = nan(size(tf_data(i).transform_abs));
        for fr=1:size(tf_data(i).transform_abs, 1)
            for k=1:n_trials
                baseline_mean = mean(tf_data(i).transform_abs(fr,baseline_range,k));
                baseline_std = std(tf_data(i).transform_abs(fr,baseline_range,k));
                tf_data(i).norm(fr,:,k) = (tf_data(i).transform_abs(fr,:,k) - baseline_mean) / (baseline_std + EPSILON); 
            end
        end

        % --- 3. STATISTICAL EVALUATION & FEATURE EXTRACTION ---
        
        % Average across frequencies for each trial
        tf_data(i).average_freq = mean(tf_data(i).norm, 1, 'omitnan');
        
        % Reshape from [1, n_times, n_trials] to [n_times, n_trials] for downstream processing
        tf_data(i).average_freq = squeeze(tf_data(i).average_freq);
        
        % Grand average across all trials (useful for visualization/plotting)
        tf_data(i).average_freqT = mean(tf_data(i).average_freq, 2, 'omitnan');

        % Pre-allocate statistical arrays
        tf_data(i).p = nan(1, length(post_baseline_range));
        tf_data(i).sign = zeros(1, length(post_baseline_range));

        for j = 1:length(post_baseline_range)
            bin_idx = post_baseline_range(j);
            
            % One-sample t-test against zero (Z-scored baseline mean is 0)
            [~, tf_data(i).p(j)] = ttest(tf_data(i).average_freq(bin_idx, :), 0);

            % Directional check: verify if the post-stimulus effect represents an activation/increase
            if mean(tf_data(i).average_freq(bin_idx, :), 'omitnan') > 0
                tf_data(i).sign(j) = 1;
            end
        end

        % Identify significant time-bins satisfying both alpha threshold and directionality
        tf_data(i).idxs = find(tf_data(i).p < pthr);
        tf_data(i).idx_sign = find(tf_data(i).sign == 1);
        tf_data(i).comparison = intersect(tf_data(i).idxs, tf_data(i).idx_sign);

        % Temporal clustering and response latency calculation
        tf_data(i).response_latency = NaN;
        if ~isempty(tf_data(i).comparison)
            % Identify discontinuities within the significant index array
            breaks = [0, find(diff(tf_data(i).comparison) ~= 1), length(tf_data(i).comparison)];
            % Cluster consecutive significant indices
            grouped = mat2cell(tf_data(i).comparison, 1, diff(breaks));
            
            for t = 1:length(grouped)
                % Validate cluster length against the minimum threshold (consecutive bin constraint)
                if length(grouped{t}) >= MIN_CLUSTER_LENGTH
                    channel_name = tf_data(i).label;
                    if ~ismember(channel_name, res_leads)
                        res_leads{end+1} = channel_name;
                    end
                    
                    % Extract onset latency from the first valid significant cluster
                    if isnan(tf_data(i).response_latency)
                        first_bin_index_relative = grouped{t}(1);
                        % Map the relative cluster index back to the absolute time vector
                        absolute_time_index = post_baseline_range(first_bin_index_relative);
                        tf_data(i).response_latency = v_times(absolute_time_index);
                        break; % Onset latency found; terminate cluster evaluation
                    end
                end
            end
        end
        
        % Workspace assignment for diagnostic tracking (legacy compatibility)
        if strcmp(condition, 'pre')
            assignin('base', 'res_leads_pre', res_leads);
        elseif strcmp(condition, 'post')
            assignin('base', 'res_leads_post', res_leads);
        end
    end

    % Post-processing: Filter out non-neural/extracranial artifact channels
    excluded_channels = ["EOG", "EKG", "DC01", "DC02", "DC03", "DC04", "DC09", "DEL SX", "DEL DX"];
    res_leads = res_leads(~ismember(string(res_leads), excluded_channels));
end