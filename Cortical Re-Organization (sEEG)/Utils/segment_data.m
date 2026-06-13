function data_3D = segment_data(data, s_event, pre_time, post_time)
    n_trials = length(s_event);
    n_samples = post_time - pre_time + 1;
    data_3D = zeros(size(data,1), n_samples, n_trials);
    
    for i = 1:n_trials
        data_3D(:,:,i) = data(:, s_event(i) + pre_time : s_event(i) + post_time);
    end
end