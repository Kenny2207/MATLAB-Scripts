%% ========================================================================
%  MAIN PIPELINE: CORTICAL REORGANIZATION ANALYSIS
%  ========================================================================
%  Target Dataset: Human sEEG (Anonymized: PAT1)
%  Context: Pre- vs. Post-Intervention Electrophysiological Modulation


%% --- GLOBAL PATH CONFIGURATION (Reproducibility Backbone) ---
% Modify this root directory path to match your local project environment setup
PROJECT_DIR = 'path/to/your/project/folder/CNR_PAT1';

%% --- DATA INGESTION & IMPORTATION ---
    EEGFILE_PRE = fullfile(PROJECT_DIR, 'Funzionali_PAT1', 'Pre', 'sub_PAT1_pre.EEG');
    EEGFILE_POST = fullfile(PROJECT_DIR, 'Funzionali_PAT1', 'Post', 'PAT1_post.EEG');
    
    [data_pre, cnames_pre, s_event_pre] = newimport_nk(EEGFILE_PRE);
    [data_post, cnames_post, s_event_post] = newimport_nk(EEGFILE_POST);
    
%% --- SIGNAL SEGMENTATION (EPOCHING) ---
    srate = 1000; % Acquisition sampling rate in Hz
    data_3D_pre = segment_data(data_pre, s_event_pre, -342, 741);
    data_3D_post = segment_data(data_post, s_event_post, -342, 741);
    
%% --- ANATOMICAL COORDINATE INTEGRATION & CHANNEL FILTERING ---
    % Load gray matter structural coordinates derived from MNI space tsv records
    tsv_path = fullfile(PROJECT_DIR, 'anatomia_PAT1', 'patientPAT1_space-MNI152NLin2009aSym.tsv');
    coord = importdata(tsv_path);
    tsv_labels = string(coord.textdata(2:end, 1));  % Extract structural labels from table rows
    thc_labels = {'U4', 'U5', 'U6', 'U7', 'U8', 'W7', 'W8', 'W9', 'Y7', 'Y8', 'Y11', 'Y12', 'Y13'};
    
    % -- Filter PRE-Intervention channels matching anatomical tsv atlas entries
    cnames_pre_str = string(cnames_pre);
    keep_idx_pre = ismember(cnames_pre_str, tsv_labels);
    data_3D_pre = data_3D_pre(keep_idx_pre, :, :);
    cnames_pre = cnames_pre(keep_idx_pre);
    
    % -- Filter POST-Intervention channels matching anatomical tsv atlas entries
    cnames_post_str = string(cnames_post);
    keep_idx_post = ismember(cnames_post_str, tsv_labels);
    data_3D_post = data_3D_post(keep_idx_post, :, :);
    cnames_post = cnames_post(keep_idx_post);
    
    % Memory optimization: clear temporary tracking variables
    clear cnames_post_str cnames_pre_str keep_idx_post keep_idx_pre

%% --- ARTIFACT IDENTIFICATION & VISUAL INSPECTION DIAGNOSTICS ---
    % Diagnostic plotting: visualize raw channel averages to isolate transients
    plot_single_channel_raw_average(data_3D_post, 26);
    plot_artifact_average(data_3D_post, 348, 10, 10);

%% --- LOCALIZED ARTIFACT REMOVAL & WAVEFORM RECONSTRUCTION ---
    disp('--- Executing PRE-intervention data correction routines ---');
    % data_3D_pre = interpolate_artifact_hardcoded(data_3D_pre);
    
    disp('--- Executing POST-intervention data correction routines ---');
    data_3D_post = interpolate_artifact_hardcoded(data_3D_post, 343:347);

%% --- TIME-FREQUENCY ANALYSIS & SPECTRAL ACTIVATION SCREENING ---
    % Define specific evaluation timescales and downsampled boundaries
    v_times_pre = linspace(47, 1042, 200);
    v_times_post = linspace(47, 1042, 200);
    baseline_range_pre = 1:60;
    post_baseline_range_pre = 61:200;
    baseline_range_post = 1:60;
    post_baseline_range_post = 61:200;

    % Execute Time-Frequency decomposition and extract responsive spectral tracks
    [tf_pre, res_leads_pre] = time_freq_analysis(data_3D_pre, cnames_pre, s_event_pre, srate, v_times_pre, 0.05 / 140, baseline_range_pre, post_baseline_range_pre, 'pre');
    [tf_post, res_leads_post] = time_freq_analysis(data_3D_post, cnames_post, s_event_post, srate, v_times_post, 0.05 / 140, baseline_range_post, post_baseline_range_post, 'post');

%% --- STATISTICAL COMPARISONS: PRE VS. POST CONDITIONS ---
    % Execute paired statistical evaluations across condition weights
    output_folder_ttest = fullfile(PROJECT_DIR, 'PAT1ttest.csv');
    ttest_result = paired_ttest_prepost(tf_pre, tf_post, res_leads_pre, res_leads_post, output_folder_ttest);

    % Map and analyze shared active contacts between conditions
    output_folder_shared = fullfile(PROJECT_DIR, 'BoxPlot_Std_Pre_Post');
    shared_results = analyze_shared_contacts(tf_pre, tf_post, cnames_pre, cnames_post, res_leads_pre, res_leads_post, output_folder_shared);

    % Compute cross-correlation of significant signal channels
    output_folder_crosscorr = fullfile(PROJECT_DIR, 'Crosscorrelazione');
    crosscorr = compute_crosscorr(tf_pre, tf_post, cnames_pre, cnames_post, shared_results, output_folder_crosscorr);

    % Track and characterize lost-response profiles during the POST condition
    output_folder_lost = fullfile(PROJECT_DIR, 'BoxPlot_Std_Lost_Contact');
    lost_response = analyze_lost_contacts(tf_pre, tf_post, cnames_pre, cnames_post, res_leads_pre, res_leads_post, output_folder_lost);

    % Standard deviation screening for baseline state evaluation
    output_folder_baseline = fullfile(PROJECT_DIR, 'BoxPlot_BASELINE_Std_Pre_Post');
    baseline_values = analyze_baseline_std(tf_pre, tf_post, cnames_pre, cnames_post, res_leads_pre, res_leads_post, output_folder_baseline);

%% --- GRAPHICAL REPRESENTATIONS & PLOTTING ROUTINES ---
    % Generate condition-specific response graphs for significant sites
    output_folder_graphs_pre = fullfile(PROJECT_DIR, 'Significant_Graphs_Pre');
    plot_significant_contacts(tf_pre, res_leads_pre, post_baseline_range_pre, 'pre', output_folder_graphs_pre);

    output_folder_graphs_post = fullfile(PROJECT_DIR, 'Significant_Graphs_Post');
    plot_significant_contacts(tf_post, res_leads_post, post_baseline_range_post, 'post', output_folder_graphs_post);

    % Generate multi-lead intersection distribution graphs
    output_folder_intersect_sig = fullfile(PROJECT_DIR, 'Significant_Graphs_pRE_Post_significant_leads');
    output_folder_intersect_all = fullfile(PROJECT_DIR, 'Significant_Graphs_pRE_Post_all_leads');
    plot_pre_post_significant(tf_pre, tf_post, res_leads_pre, res_leads_post, shared_results, ttest_result, output_folder_intersect_sig, output_folder_intersect_all);
 
    % Single channel targeted plotting routine (Example lead: Q18)
    canale_desiderato = "Q18"; 
    limite_Y_min = -0.2;          
    limite_Y_max = 0.5;           
    output_folder_single_lost = fullfile(PROJECT_DIR, 'SingleChannel_LOST');
    plot_single_channel_significant(tf_pre, tf_post, ttest_result, canale_desiderato, limite_Y_min, limite_Y_max, output_folder_single_lost, true);

    % Automated axis-scaling plotting routine (Example lead: J14)
    output_dir_single_auto = fullfile(PROJECT_DIR, 'Singolo_canale');
    canale_da_plot = 'J14'; 
    plot_sem = true; 
    plot_single_channel_auto_ylim(tf_pre, tf_post, ttest_result, output_dir_single_auto, canale_da_plot, plot_sem);

%% --- INTRACRANIAL EVOKED RESPONSES (iERP) PIPELINE ---
    v_times = linspace(-300, 699, 200);  
    
    % Execute iERP processing for the PRE-intervention condition
    [erp_data_pre, res_leads_ierp_pre, v_times_ierp] = compute_iERP(data_3D_pre, cnames_pre, 1000, 0.1, 35, 3, 0.05/140, 1:60, 61:200, false);
    output_folder_ierp_pre = fullfile(PROJECT_DIR, 'Significant_Graphs_iERP_pre');
    plot_significant_iERP(erp_data_pre, res_leads_ierp_pre, v_times_ierp, 'pre', output_folder_ierp_pre);

    % Execute iERP processing for the POST-intervention condition
    [erp_data_post, res_leads_ierp_post, v_times_ierp_post] = compute_iERP(data_3D_post, cnames_post, srate, 0.1, 35, 3, 0.05/140, 1:60, 61:200, false);
    output_folder_ierp_post = fullfile(PROJECT_DIR, 'Significant_Graphs_iERP_post');
    plot_significant_iERP(erp_data_post, res_leads_ierp_post, v_times_ierp_post, 'post', output_folder_ierp_post);

    % Extract and analyze shared tracking leads across evoked datasets
    output_folder_ierp_shared = fullfile(PROJECT_DIR, 'Significant_Graphs_iERP_shared');
    shared_ierp_results = analyze_shared_ierp_results(erp_data_pre, erp_data_post, cnames_pre, cnames_post, res_leads_ierp_pre, res_leads_ierp_post, output_folder_ierp_shared);

    % Execute paired t-tests on iERP timeseries records
    ttest_ierp_results = paired_ttest_ierp_prepost(erp_data_pre, erp_data_post, res_leads_ierp_pre, res_leads_ierp_post);

    % Plot intersection maps for combined evoke responses
    output_folder_ierp_intersect = fullfile(PROJECT_DIR, 'Significant_Graphs_iERP_pRE_Post');
    plot_pre_post_ierp_significant(erp_data_pre, erp_data_post, res_leads_ierp_pre, res_leads_ierp_post, shared_ierp_results, ttest_ierp_results, output_folder_ierp_intersect);

%% --- METADATA EXPORTATION & SHEET GENERATION ---
    output_report_path = fullfile(PROJECT_DIR, 'CNR_PAT1_FULL_REPORT.xlsx');
    export_full_report('PAT1', thc_labels, res_leads_pre, res_leads_post, ...
                       shared_results, crosscorr, lost_response, ...
                       res_leads_ierp_pre, res_leads_ierp_post, ...
                       shared_ierp_results, ttest_ierp_results, ...
                       output_report_path);
%% ========================================================================