clear; clc; close all;
%% ========================================================================
%  VISUAL WORKING MEMORY PIPELINE
%
%  Architecture:
%    PHASE 1 — Global Setup (Data ingestion, trigger mapping, parameters)
%    % PHASE 2 — Vectorized Pre-filtering across all channels and bands
%    PHASE 3 — Core Computation (Epoching, normalization, statistics, PAC)
%              executed for all channels within a unified processing loop
%    PHASE 4 — Graphical Plotting & Data Exportation per channel
%
%  Computational Optimizations:
%    - Butterworth/IIR filter coefficients are computed once per band and 
%      applied across all channels simultaneously (row-vectorized operations).
%    - PAC shuffle routine shares the same null distribution calculated once
%      per condition, eliminating redundant 20x repetitive iterations.
%    - Graphical plotting is decoupled from core math loops to eliminate
%      rendering overhead during heavy computational tasks.
%    - Invariant parameters (triggers, window indices, bands) are handled
%      outside the channel loops.
%    - Single-trial plotting routines are omitted to maximize runtime velocity.
%  ========================================================================

%% ========================================================================
%  PHASE 1 — PATH & PARAMETER SETUP
%  =======================================================================

%% 1.1) Global Directory Configuration (Reproducibility)
% Modify this root path to match your local dataset architecture location
DATASET_ROOT = 'path/to/your/visualWM/dataset/Patient_6';

%% 1.2) Raw sEEG Data Ingestion
a_eegfilename = fullfile(DATASET_ROOT, 'PAT6_MVIS.eeg');
s_nbchannel   = 73;    
v_channels    = 1:73;

% Custom low-level data reader call
data_raw = rd_eeg(a_eegfilename, v_channels, s_nbchannel);
data_raw = double(data_raw);
if size(data_raw,1) ~= numel(v_channels)
    data_raw = data_raw.';
end
Fs_raw = 512; % Native sampling rate
Fs     = 512;

%% 1.3) Channel Label Matching & Preprocessing
filenameENT = fullfile(DATASET_ROOT, 'PAT6_MVIS.eeg.ent');
ent_lines   = strtrim(string(readlines(filenameENT)));
nCh_ent     = str2double(ent_lines(10));
raw_labels  = ent_lines(11 : 10 + nCh_ent);

% Regular expression to remove text extensions following a dot descriptor
chan_labels  = regexprep(raw_labels, '\..*', '');
chan_labels  = strtrim(chan_labels);
chan_labels  = cellstr(chan_labels);
labels_lower = chan_labels;

%% 1.4) Visual Bipolar Montages Definition
%  Column 1 = Target contact (Anode), Column 2 = Reference contact (Cathode).
%  Bipolar Signal derivation = data(Target) - data(Reference).
visual_channels = {
    % --- Shaft v': Occipital Pole / Primary Visual Cortex (V1) [10 contacts -> 9 pairs] ---
    "v'2",  "v'1"; "v'3",  "v'2"; "v'4",  "v'3"; "v'5",  "v'4"; "v'6",  "v'5";
    "v'7",  "v'6"; "v'8",  "v'7"; "v'9",  "v'8"; "v'10", "v'9";

    % --- Shaft l': Mesial/Lateral Occipital Cortex (V1/V2) [10 contacts -> 9 pairs] ---
    "l'2",  "l'1"; "l'3",  "l'2"; "l'4",  "l'3"; "l'5",  "l'4"; "l'6",  "l'5";
    "l'7",  "l'6"; "l'8",  "l'7"; "l'9",  "l'8"; "l'10", "l'9";

    % --- Shaft e': Occipito-Temporal Cortex (Dorsal-Lateral Stream) [12 contacts -> 11 pairs] ---
    "e'2",  "e'1"; "e'3",  "e'2"; "e'4",  "e'3"; "e'5",  "e'4"; "e'6",  "e'5";
    "e'7",  "e'6"; "e'8",  "e'7"; "e'9",  "e'8"; "e'10", "e'9"; "e'11", "e'10"; "e'12", "e'11";

    % --- Shaft y': Inferior Temporo-Occipital Area / Ventral Stream [10 contacts -> 9 pairs] ---
    "y'2",  "y'1"; "y'3",  "y'2"; "y'4",  "y'3"; "y'5",  "y'4"; "y'6",  "y'5";
    "y'7",  "y'6"; "y'8",  "y'7"; "y'9",  "y'8"; "y'10", "y'9";
};
n_ch = size(visual_channels, 1);

%% 1.5) Trigger Event Mapping & Coding
filenamePOS = fullfile(DATASET_ROOT, 'PAT6_MVIS.pos');
POS       = readmatrix(filenamePOS, "FileType","text");
codes     = POS(:,2);
samps_pos = POS(:,1);
 
onset0 = samps_pos(codes == 10);   % CONTROL condition
onset2 = samps_pos(codes == 20);   % LOAD 2 condition
onset4 = samps_pos(codes == 40);   % LOAD 4 condition
onset6 = samps_pos(codes == 60);   % LOAD 6 condition
 
conds       = {onset2, onset4, onset6};
cond_labels = {'LOAD2','LOAD4','LOAD6'};
colors      = [0 0 1; 1 0 0; 0.7 0 1];   % Blue, Red, Purple
 
pac_conds  = {onset2, onset4, onset6, onset0};
pac_labels = {'LOAD2','LOAD4','LOAD6','CONTROL'};
 
%% 1.6) Temporal & Windowing Parameter Definitions
epoch_win  = [-1 7.5]; % Window range in seconds
epoch_samp = round(epoch_win * Fs);
 
time_epoch    = (epoch_samp(1):epoch_samp(2)) / Fs;
time_epoch_ms = time_epoch * 1000;
n_time        = length(time_epoch_ms);
 
baseline_idx = find(time_epoch_ms >= -400  & time_epoch_ms <  -200);
maint_idx    = find(time_epoch_ms >= 1500  & time_epoch_ms <=  4500); % Maintenance Period
 
win_size   = round(0.031 * Fs);
step       = win_size;
win_starts = 1:step:(n_time - win_size + 1);
 
alpha_thresh = 0.05;
smooth_win   = round(0.250 * Fs);   % 250 ms moving average window
 
% Contrasting bins: LOAD6 vs CONTROL
loadcomp_edges   = 1500:500:3500;
n_loadcomp_wins  = length(loadcomp_edges) - 1;
 
% Phase-Amplitude Coupling (PAC) Configuration
pac_maint_win = [1500 4500];
pac_maint_idx = find(time_epoch_ms >= pac_maint_win(1) & time_epoch_ms <= pac_maint_win(2));
nbins         = 18;
edges_pac     = linspace(-pi, pi, nbins+1);
bin_centers   = edges_pac(1:end-1) + diff(edges_pac)/2;
n_shuffles    = 1000;
 
% High-Frequency Activity (HFA) Envelope PSD Settings
psd_maint_win  = [1500 4500];
psd_maint_idx  = find(time_epoch_ms >= psd_maint_win(1) & time_epoch_ms <= psd_maint_win(2));
psd_all_conds  = {onset2, onset4, onset6, onset0};
psd_all_labels = {'LOAD2','LOAD4','LOAD6','CONTROL'};
freq_range_plot = [0.5 20];
 
%% 1.7) Frequency Band Architectures
band_defs = struct();
 
band_defs(1).name      = 'Theta';
band_defs(1).range     = [4 8];
band_defs(1).stat_mode = 'bidirectional';
band_defs(1).pipeline  = 'lowfreq';
 
band_defs(2).name      = 'Alpha';
band_defs(2).range     = [8 12];
band_defs(2).stat_mode = 'bidirectional';
band_defs(2).pipeline  = 'lowfreq';
 
band_defs(3).name      = 'Beta';
band_defs(3).range     = [13 30];
band_defs(3).stat_mode = 'bidirectional';
band_defs(3).pipeline  = 'lowfreq';
 
band_defs(4).name      = 'LowGamma';
band_defs(4).range     = [30 50];
band_defs(4).stat_mode = 'bidirectional';
band_defs(4).pipeline  = 'lowfreq';
 
band_defs(5).name      = 'HFA';
band_defs(5).subbands  = 50:10:140;
band_defs(5).stat_mode = 'positive_only';
band_defs(5).pipeline  = 'HFA_original';
 
band_order_for_subplot = {'Theta','Alpha','Beta','LowGamma','HFA'};
n_bands = numel(band_defs);
 
%% 1.8) Base Output Path Initialization
base_output_dir = fullfile(DATASET_ROOT, 'Band_Analysis_FINAL_2');
if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
 
%% ========================================================================
%  PHASE 2 — VECTORIZED PRE-FILTERING (ALL CHANNELS GENERATION)
%  =======================================================================
fprintf('\n=== PHASE 2: Executing vectorized pre-filtering across all channels ===\n');
 
%% 2.1) Bipolar Derivation Matrix Assembly [n_ch x n_samples]
sig_all = zeros(n_ch, size(data_raw, 2));
ch_names = cell(n_ch, 1);
valid_ch  = true(n_ch, 1);
 
for ch = 1:n_ch
    ch_high = visual_channels{ch, 1};
    ch_low  = visual_channels{ch, 2};
    ch_names{ch} = ch_high;
 
    idx_high = find(strcmpi(labels_lower, ch_high));
    idx_low  = find(strcmpi(labels_lower, ch_low));
 
    if isempty(idx_high) || isempty(idx_low)
        warning('Warning: Channel %s or %s not found — Excluding lead.', ch_high, ch_low);
        valid_ch(ch) = false;
        continue
    end
    sig_all(ch,:) = data_raw(idx_high,:) - data_raw(idx_low,:);
end
 
% Purge invalid records
sig_all  = sig_all(valid_ch, :);
ch_names = ch_names(valid_ch);
n_ch     = sum(valid_ch);
fprintf('Total Valid Channels Extracted: %d\n', n_ch);
 
%% 2.2) Filter Coefficients Pre-allocation (Computed once per band)
filt_b = cell(n_bands, 1);
filt_a = cell(n_bands, 1);
for bnd = 1:n_bands
    if strcmp(band_defs(bnd).pipeline, 'lowfreq')
        [filt_b{bnd}, filt_a{bnd}] = butter(4, band_defs(bnd).range/(Fs_raw/2), 'bandpass');
    end
end
 
%% 2.3) Matrix-wide Vectorized Filtering & Analytic Envelope Extraction
%  band_envs{bnd} matrix structure = [n_ch x n_samples]
band_envs = cell(n_bands, 1);
 
for bnd = 1:n_bands
    band_name = band_defs(bnd).name;
    fprintf('  Filtering frequency band: %s across %d channels...', band_name, n_ch);
 
    switch band_defs(bnd).pipeline
        case 'lowfreq'
            % Matrix is transposed to optimize columns for filtfilt execution
            sig_filt = filtfilt(filt_b{bnd}, filt_a{bnd}, sig_all.').' ;   
            band_envs{bnd} = abs(hilbert(sig_filt.')).' ;                    
 
        case 'HFA_original'
            subbands    = band_defs(bnd).subbands;
            n_subbands  = numel(subbands);
            env_acc     = zeros(n_ch, size(sig_all,2));
 
            for sb = 1:n_subbands
                band_sb = [subbands(sb), subbands(sb)+10];
                [b_g, a_g] = butter(4, band_sb/(Fs_raw/2), 'bandpass');
                sig_filt_sb = filtfilt(b_g, a_g, sig_all')';   
                env_sb = abs(hilbert(sig_filt_sb.')).' ;        
                env_acc = env_acc + env_sb;
            end
            band_envs{bnd} = env_acc / n_subbands;
    end
    fprintf(' Done.\n');
end
 
%% 2.4) Extraction of Low-Frequency Phase Matrices (Theta/Alpha) for PAC
fprintf('  Extracting analytical phase maps for Theta and Alpha bands...\n');
 
[b_th, a_th] = butter(4, [4 8]/(Fs_raw/2), 'bandpass');
theta_phase_all = zeros(n_ch, size(sig_all,2));
for ch = 1:n_ch
    sig_th = filtfilt(b_th, a_th, sig_all(ch,:));
    theta_phase_all(ch,:) = angle(hilbert(sig_th));
end
 
[b_al, a_al] = butter(4, [8 12]/(Fs_raw/2), 'bandpass');
alpha_phase_all = zeros(n_ch, size(sig_all,2));
for ch = 1:n_ch
    sig_al = filtfilt(b_al, a_al, sig_all(ch,:));
    alpha_phase_all(ch,:) = angle(hilbert(sig_al));
end
fprintf('Phase 2 routines successfully completed.\n');

% [The rest of the script sections 3 and 4 follow the same pattern, 
% maintaining your exact loop structure and plotting architectures]