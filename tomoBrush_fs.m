%% Script to load 2 WAV files, find peaks, extract segments, and perform pspectrum Analysis
% With additional feature extraction in 4000-10000 Hz range

% Clear workspace and close all figures
clear all;
close all;

% Parameters for peak detection and segment extraction
target_peaks = 5;              % Number of peaks to detect
smooth_window_size = 0.05;     % Window size for smoothing in seconds
min_peak_distance_slow = 0.5;  % Minimum distance between peaks for slow file
min_peak_distance_fast = 0.2;  % Minimum distance between peaks for fast file
min_peak_prominence = 0.1;     % Minimum peak prominence as fraction of max
segment_duration = 0.2;        % Duration of segment to extract around each peak (in seconds)

% Feature extraction parameters
freq_range_min = 4000;         % Minimum frequency for feature extraction (Hz)
freq_range_max = 10000;        % Maximum frequency for feature extraction (Hz)
num_features = 10;             % Number of features to extract from the specified range

% File names
file1 = 'TopMidBot_Slow5.wav';
file2 = 'TopMidBot_Fast5.wav';
file_pairs = {file1(1:end-4), file2(1:end-4)};  % Remove .wav extension

% Process each file
slow_segments = [];
fast_segments = [];
feature_data = struct();  % Structure to store extracted features

for i = 1:length(file_pairs)
    file_prefix = file_pairs{i};
    fprintf('\nProcessing file: %s.wav\n', file_prefix);
    
    %% Load Audio
    [audio, fs] = audioread([file_prefix, '.wav']);
    if size(audio,2) > 1
        audio = mean(audio,2);  % Convert to mono if stereo
    end
    audio = audio / max(abs(audio)); % Normalize audio
    t = (0:length(audio)-1) / fs;
    
    %% Envelope Calculation and Smoothing
    env = abs(audio);  % Simple envelope - absolute value
    smooth_window_samples = round(smooth_window_size * fs);
    if mod(smooth_window_samples,2) == 0
        smooth_window_samples = smooth_window_samples + 1;  % ensure odd window size
    end
    smooth_kernel = ones(smooth_window_samples,1) / smooth_window_samples;
    smooth_env = conv(env, smooth_kernel, 'same');
    
    %% Peak Detection
    % Choose min peak distance (in seconds) based on file type:
    if contains(file_prefix, 'Slow')
        min_peak_dist_sec = min_peak_distance_slow;
    else
        min_peak_dist_sec = min_peak_distance_fast;
    end
    min_peak_dist_samples = round(min_peak_dist_sec * fs);
    
    % Detect peaks using findpeaks with minimum distance and prominence
    [peak_values, peak_indices] = findpeaks(smooth_env, 'MinPeakDistance', min_peak_dist_samples, ...
        'MinPeakProminence', min_peak_prominence * max(smooth_env));
    
    if isempty(peak_indices)
        fprintf('Warning: No peaks detected in %s.wav\n', file_prefix);
        continue;
    end
    
    % Sort detected peaks in descending order (largest first)
    [~, order] = sort(peak_values, 'descend');
    num_detected = length(order);
    
    if num_detected < target_peaks
        fprintf('Warning: Only %d peaks detected in %s, expected %d. Using all detected peaks.\n', ...
            num_detected, file_prefix, target_peaks);
        selected_order = order;
    else
        selected_order = order(1:target_peaks);
    end
    
    % Sort the selected peak indices in chronological order
    selected_indices = sort(peak_indices(selected_order));
    
    % Store peak times and values
    peak_times = t(selected_indices);
    selected_values = smooth_env(selected_indices);
    
    %% Segment Extraction around Each Peak
    half_seg_samples = round((segment_duration/2) * fs);
    num_peaks = length(selected_indices);
    segments = zeros(num_peaks, 2*half_seg_samples + 1);  % Preallocate for each segment
    
    for k = 1:num_peaks
        center = selected_indices(k);
        idx_range = (center-half_seg_samples):(center+half_seg_samples);
        % Ensure indices remain within bounds
        idx_range = idx_range(idx_range >= 1 & idx_range <= length(audio));
        seg = audio(idx_range);
        % Zero-pad if segment is shorter than expected
        if length(seg) < (2*half_seg_samples+1)
            pad_length = (2*half_seg_samples+1) - length(seg);
            seg = [seg; zeros(pad_length, 1)];
        end
        segments(k,:) = seg;
    end
    
    %% pspectrum Analysis of Segments
    % For each segment, compute its power spectrum using pspectrum
    [p_temp, f_ps] = pspectrum(segments(1,:), fs);
    num_freq = length(p_temp);
    pspec_results = zeros(num_peaks, num_freq);  % Preallocate
    
    figure(i+2000)
    hold on;
    for k = 1:num_peaks
        [p_segment, f_ps] = pspectrum(segments(k,:), fs);
        pspec_results(k,:) = p_segment;
        plot(f_ps, p_segment, 'Color', [0.7, 0.7, 0.7]);  % Plot individual pspectra in light gray
    end
    
    avg_ps = mean(pspec_results, 1);  % Average over segments
    env_avg_ps = envelope(avg_ps);    % Compute envelope of the mean spectrum
    
    plot(f_ps, env_avg_ps, 'r', 'LineWidth', 2)  % Overlay in red
    title(sprintf('Envelope + pspectrum for Segment Group %d', i))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    legend('Individual pspectra', 'Envelope of Mean pspectrum')
    hold off;
    
    %% Feature Extraction from 4000-10000 Hz range
    % Find indices for our frequency range of interest
    freq_indices = find(f_ps >= freq_range_min & f_ps <= freq_range_max);
    
    if isempty(freq_indices)
        fprintf('Warning: No frequencies found in the specified range (%.0f-%.0f Hz)\n', ...
            freq_range_min, freq_range_max);
        continue;
    end
    
    % Extract frequency and power values in our range of interest
    selected_freqs = f_ps(freq_indices);
    selected_powers = avg_ps(freq_indices);
    
    % Create a new figure to display the selected frequency range
    figure(i+3000)
    plot(selected_freqs, selected_powers, 'b', 'LineWidth', 1.5);
    title(sprintf('Selected Frequency Range (%.0f-%.0f Hz) for %s', ...
        freq_range_min, freq_range_max, file_prefix));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    
    % Feature 1: Total energy in the selected frequency band
    feature_data.(file_prefix).total_energy = sum(selected_powers);
    
    % Feature 2: Mean energy in the selected frequency band
    feature_data.(file_prefix).mean_energy = mean(selected_powers);
    
    % Feature 3: Spectral centroid within the selected band
    feature_data.(file_prefix).centroid = sum(selected_freqs .* selected_powers) / sum(selected_powers);
    
    % Feature 4: Spectral flatness (ratio of geometric mean to arithmetic mean)
    % Add small value to avoid log(0) problems
    epsilon = 1e-10;
    log_powers = log(selected_powers + epsilon);
    feature_data.(file_prefix).flatness = exp(mean(log_powers)) / mean(selected_powers);
    
    % Feature 5: Spectral variance
    feature_data.(file_prefix).variance = var(selected_powers);
    
    % Feature 6: Find peak frequencies in selected range
    [peak_vals, peak_locs] = findpeaks(selected_powers, 'SortStr', 'descend');
    
    if length(peak_vals) >= 3
        % Store the top 3 peak frequencies
        top_peak_freqs = selected_freqs(peak_locs(1:3));
        feature_data.(file_prefix).peak_freqs = top_peak_freqs;
        
        % Feature 7: Peak-to-average ratio
        feature_data.(file_prefix).peak_avg_ratio = peak_vals(1) / mean(selected_powers);
    else
        % If fewer than 3 peaks, store what we have
        top_peak_freqs = selected_freqs(peak_locs(1:min(3,length(peak_locs))));
        feature_data.(file_prefix).peak_freqs = top_peak_freqs;
        
        if ~isempty(peak_vals)
            feature_data.(file_prefix).peak_avg_ratio = peak_vals(1) / mean(selected_powers);
        else
            feature_data.(file_prefix).peak_avg_ratio = 0;
        end
    end
    
    % Feature 8: Energy distribution (divide band into sub-bands)
    num_subbands = 5;
    subband_indices = round(linspace(1, length(selected_freqs), num_subbands+1));
    subband_energies = zeros(1, num_subbands);
    
    for sb = 1:num_subbands
        idx_range = subband_indices(sb):subband_indices(sb+1);
        if length(idx_range) > 1  % Ensure we have at least 2 points
            subband_energies(sb) = sum(selected_powers(idx_range));
        end
    end
    
    feature_data.(file_prefix).subband_energies = subband_energies;
    feature_data.(file_prefix).subband_energy_ratio = subband_energies / sum(subband_energies);
    
    % Feature 9: Energy slope (linear regression on spectrum)
    x = selected_freqs;
    y = selected_powers;
    p = polyfit(x, y, 1);
    feature_data.(file_prefix).energy_slope = p(1);
    
    % Feature 10: Energy roll-off point (frequency below which 85% of energy is contained)
    cumulative_energy = cumsum(selected_powers);
    roll_off_point = find(cumulative_energy >= 0.85 * cumulative_energy(end), 1, 'first');
    if ~isempty(roll_off_point)
        feature_data.(file_prefix).roll_off_freq = selected_freqs(roll_off_point);
    else
        feature_data.(file_prefix).roll_off_freq = NaN;
    end
    
    %% Plot results from earlier steps
    figure('Position',[100,100,1200,900]);
    
    % Plot 1: Waveform with detected peaks
    subplot(3,1,1);
    plot(t, audio); hold on;
    plot(peak_times, audio(selected_indices), 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
    title(sprintf('Waveform with Detected Peaks: %s', file_prefix));
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
    
    % Plot 2: Example segments (first 5)
    subplot(3,1,2);
    segment_time = (-half_seg_samples:half_seg_samples) / fs;
    hold on;
    for k = 1:min(5, num_peaks)
        plot(segment_time, segments(k,:) + (k-1)*1.5);  % Offset for visibility
    end
    title('First 5 Extracted Segments (Offset for Visibility)');
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
    
    % Plot 3: Selected Frequency Range with Features
    subplot(3,1,3);
    plot(selected_freqs, selected_powers, 'b', 'LineWidth', 1.5); hold on;
    
    % Mark peak frequencies
    if ~isempty(feature_data.(file_prefix).peak_freqs)
        for p_idx = 1:length(feature_data.(file_prefix).peak_freqs)
            peak_freq = feature_data.(file_prefix).peak_freqs(p_idx);
            [~, idx] = min(abs(selected_freqs - peak_freq));
            plot(peak_freq, selected_powers(idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            text(peak_freq, selected_powers(idx)*1.1, sprintf('%.0f Hz', peak_freq), ...
                'FontSize', 8, 'HorizontalAlignment', 'center');
        end
    end
    
    % Mark spectral centroid
    centroid = feature_data.(file_prefix).centroid;
    % Create a vertical line at the centroid frequency
    % FIX: Add checks to prevent plotting errors
    if ~isempty(selected_powers) && isscalar(centroid) && ~isnan(centroid)
        max_power = max(selected_powers);
        if ~isnan(max_power) && max_power > 0
            plot([centroid, centroid], [0, max_power*0.8], 'g--', 'LineWidth', 1.5);
            text(centroid, max_power*0.4, sprintf('Centroid: %.0f Hz', centroid), ...
                'FontSize', 8, 'HorizontalAlignment', 'right');
        end
    end
    
    % Mark roll-off point
    roll_off = feature_data.(file_prefix).roll_off_freq;
    if ~isnan(roll_off)
        % Create a vertical line at the roll-off frequency
        % Using the maximum power for consistent scaling
        max_power = max(selected_powers);
        plot([roll_off, roll_off], [0, max_power*0.8], 'm--', 'LineWidth', 1.5);
        text(roll_off, max_power*0.3, sprintf('Roll-off: %.0f Hz', roll_off), ...
            'FontSize', 8, 'HorizontalAlignment', 'left');
    end
    
    title(sprintf('Selected Frequency Range (%.0f-%.0f Hz) with Features', ...
        freq_range_min, freq_range_max));
    xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;
    
    %% Store results for later comparison
    if contains(file_prefix, 'Slow')
        slow_peak_indices = selected_indices;
        slow_peak_times = peak_times;
        slow_peak_values = selected_values;
        slow_intervals = diff(peak_times);
        slow_segments = segments;
        slow_pspec_results = pspec_results;
        slow_avg_ps = avg_ps;        % Average pspectrum for slow segments
        slow_ps_freq = f_ps;         % Frequency axis from pspectrum
        slow_fs = fs;
        slow_selected_freqs = selected_freqs;
        slow_selected_powers = selected_powers;
    else
        fast_peak_indices = selected_indices;
        fast_peak_times = peak_times;
        fast_peak_values = selected_values;
        fast_intervals = diff(peak_times);
        fast_segments = segments;
        fast_pspec_results = pspec_results;
        fast_avg_ps = avg_ps;        % Average pspectrum for fast segments
        fast_ps_freq = f_ps;         % Frequency axis from pspectrum
        fast_fs = fs;
        fast_selected_freqs = selected_freqs;
        fast_selected_powers = selected_powers;
    end
    
    %% Display extracted features
    fprintf('\nExtracted Features for %s in %.0f-%.0f Hz range:\n', file_prefix, freq_range_min, freq_range_max);
    fprintf('  Total energy: %.4e\n', feature_data.(file_prefix).total_energy);
    fprintf('  Mean energy: %.4e\n', feature_data.(file_prefix).mean_energy);
    fprintf('  Spectral centroid: %.2f Hz\n', feature_data.(file_prefix).centroid);
    fprintf('  Spectral flatness: %.4f\n', feature_data.(file_prefix).flatness);
    fprintf('  Spectral variance: %.4e\n', feature_data.(file_prefix).variance);
    fprintf('  Peak-to-average ratio: %.4f\n', feature_data.(file_prefix).peak_avg_ratio);
    fprintf('  Energy slope: %.4e\n', feature_data.(file_prefix).energy_slope);
    fprintf('  Roll-off frequency: %.2f Hz\n', feature_data.(file_prefix).roll_off_freq);
    
    % Print peak frequencies
    if ~isempty(feature_data.(file_prefix).peak_freqs)
        fprintf('  Peak frequencies: ');
        for p_idx = 1:length(feature_data.(file_prefix).peak_freqs)
            fprintf('%.2f Hz ', feature_data.(file_prefix).peak_freqs(p_idx));
        end
        fprintf('\n');
    end
    
    % Print subband energy distribution
    fprintf('  Subband energy ratios: ');
    for sb = 1:num_subbands
        fprintf('Band %d: %.4f  ', sb, feature_data.(file_prefix).subband_energy_ratio(sb));
    end
    fprintf('\n');
end

%% Check that segments have been extracted for both recordings
if isempty(slow_segments) || isempty(fast_segments)
    error('Not enough segments extracted. Adjust your peak detection parameters.');
end

%% Final comparison plot of selected frequency ranges
if exist('slow_selected_freqs', 'var') && exist('fast_selected_freqs', 'var')
    figure('Position',[100,100,1200,600]);
    
    subplot(2,1,1);
    plot(slow_selected_freqs, slow_selected_powers, 'b', 'LineWidth',1.5);
    title(sprintf('Selected Frequency Range (%.0f-%.0f Hz) - Slow Recording', ...
        freq_range_min, freq_range_max));
    xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;
    
    subplot(2,1,2);
    plot(fast_selected_freqs, fast_selected_powers, 'r', 'LineWidth',1.5);
    title(sprintf('Selected Frequency Range (%.0f-%.0f Hz) - Fast Recording', ...
        freq_range_min, freq_range_max));
    xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;
end

%% Feature comparison between slow and fast recordings
figure('Position',[100,100,1200,800]);

% Common feature names for plotting
feature_names = {'Total Energy', 'Mean Energy', 'Centroid', 'Flatness', ...
    'Variance', 'Peak/Avg Ratio', 'Energy Slope', 'Roll-off Freq'};

% Extract common features for comparison
% Check and handle features that might not be scalar
feature_names_full = {'total_energy', 'mean_energy', 'centroid', 'flatness', 'variance', 'peak_avg_ratio', 'energy_slope', 'roll_off_freq'};
slow_features = zeros(1, length(feature_names_full));
fast_features = zeros(1, length(feature_names_full));

for i = 1:length(feature_names_full)
    % Get slow feature
    if isfield(feature_data, file_pairs{1}) && isfield(feature_data.(file_pairs{1}), feature_names_full{i})
        feat = feature_data.(file_pairs{1}).(feature_names_full{i});
        if ~isempty(feat) && isscalar(feat) && ~isnan(feat)
            slow_features(i) = feat;
        else
            slow_features(i) = 0;  % Default value if feature is invalid
            fprintf('Warning: Feature %s for slow recording is not a valid scalar\n', feature_names_full{i});
        end
    else
        slow_features(i) = 0;
        fprintf('Warning: Feature %s for slow recording does not exist\n', feature_names_full{i});
    end
    
    % Get fast feature
    if isfield(feature_data, file_pairs{2}) && isfield(feature_data.(file_pairs{2}), feature_names_full{i})
        feat = feature_data.(file_pairs{2}).(feature_names_full{i});
        if ~isempty(feat) && isscalar(feat) && ~isnan(feat)
            fast_features(i) = feat;
        else
            fast_features(i) = 0;  % Default value if feature is invalid
            fprintf('Warning: Feature %s for fast recording is not a valid scalar\n', feature_names_full{i});
        end
    else
        fast_features(i) = 0;
        fprintf('Warning: Feature %s for fast recording does not exist\n', feature_names_full{i});
    end
end

% Normalize features for better visualization (each feature between 0-1)
combined_features = [slow_features; fast_features];
min_vals = min(combined_features, [], 1);
max_vals = max(combined_features, [], 1);
range_vals = max_vals - min_vals;

% Avoid division by zero
range_vals(range_vals == 0) = 1;

% Normalize
norm_slow_features = (slow_features - min_vals) ./ range_vals;
norm_fast_features = (fast_features - min_vals) ./ range_vals;

% Create subplot for each feature
num_plot_features = length(feature_names);
for i = 1:num_plot_features
    subplot(2, 4, i);
    bar([norm_slow_features(i), norm_fast_features(i)]);
    title(feature_names{i});
    set(gca, 'XTickLabel', {'Slow', 'Fast'});
    ylabel('Normalized Value');
    grid on;
    
    % Add text showing the actual values
    text(1, norm_slow_features(i) + 0.05, sprintf('%.2e', slow_features(i)), ...
        'HorizontalAlignment', 'center');
    text(2, norm_fast_features(i) + 0.05, sprintf('%.2e', fast_features(i)), ...
        'HorizontalAlignment', 'center');
end

sgtitle(sprintf('Feature Comparison in %.0f-%.0f Hz Range', freq_range_min, freq_range_max));

%% Calculate feature differences and ratios
fprintf('\n\nFeature Comparison between Slow and Fast Recordings:\n');
fprintf('%-20s %-15s %-15s %-15s %-15s\n', 'Feature', 'Slow', 'Fast', 'Difference', 'Ratio (Fast/Slow)');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:length(feature_names)
    slow_val = slow_features(i);
    fast_val = fast_features(i);
    diff_val = fast_val - slow_val;
    
    % Handle division by zero for ratio
    if slow_val ~= 0
        ratio_val = fast_val / slow_val;
    else
        ratio_val = NaN;
    end
    
    fprintf('%-20s %-15.4e %-15.4e %-15.4e %-15.4f\n', ...
        feature_names{i}, slow_val, fast_val, diff_val, ratio_val);
end

%% Print final statistics on intervals
if ~isempty(slow_intervals) && ~isempty(fast_intervals)
    fprintf('\nComparison of average intervals:\n');
    fprintf('Slow recording: %.4f seconds\n', mean(slow_intervals));
    fprintf('Fast recording: %.4f seconds\n', mean(fast_intervals));
    fprintf('Ratio (Slow/Fast): %.2f\n', mean(slow_intervals)/mean(fast_intervals));
end

%% Classification example: Can we distinguish between slow and fast segments?
% Create feature matrix (each row is a segment, each column is a feature)
fprintf('\n\nClassification Example: Distinguishing Slow vs Fast Segments\n');

% For simplicity, we'll use just 3 distinctive features
selected_features = [3, 4, 8];  % Centroid, Flatness, Roll-off
feature_names_selected = feature_names(selected_features);

% Extract segment-level features for each selected peak
slow_segment_features = zeros(size(slow_segments, 1), length(selected_features));
fast_segment_features = zeros(size(fast_segments, 1), length(selected_features));

% Process each slow segment
for k = 1:size(slow_segments, 1)
    segment = slow_segments(k, :);
    [p_segment, f_segment] = pspectrum(segment, slow_fs);
    
    % Find indices for our frequency range
    freq_idx = find(f_segment >= freq_range_min & f_segment <= freq_range_max);
    
    if ~isempty(freq_idx)
        sel_freq = f_segment(freq_idx);
        sel_power = p_segment(freq_idx);
        
        % Centroid
        slow_segment_features(k, 1) = sum(sel_freq .* sel_power) / sum(sel_power);
        
        % Flatness
        epsilon = 1e-10;
        log_powers = log(sel_power + epsilon);
        slow_segment_features(k, 2) = exp(mean(log_powers)) / mean(sel_power);
        
        % Roll-off
        cumulative_energy = cumsum(sel_power);
        roll_off_point = find(cumulative_energy >= 0.85 * cumulative_energy(end), 1, 'first');
        if ~isempty(roll_off_point)
            slow_segment_features(k, 3) = sel_freq(roll_off_point);
        else
            slow_segment_features(k, 3) = NaN;
        end
    end
end

% Process each fast segment
for k = 1:size(fast_segments, 1)
    segment = fast_segments(k, :);
    [p_segment, f_segment] = pspectrum(segment, fast_fs);
    
    % Find indices for our frequency range
    freq_idx = find(f_segment >= freq_range_min & f_segment <= freq_range_max);
    
    if ~isempty(freq_idx)
        sel_freq = f_segment(freq_idx);
        sel_power = p_segment(freq_idx);
        
        % Centroid
        fast_segment_features(k, 1) = sum(sel_freq .* sel_power) / sum(sel_power);
        
        % Flatness
        epsilon = 1e-10;
        log_powers = log(sel_power + epsilon);
        fast_segment_features(k, 2) = exp(mean(log_powers)) / mean(sel_power);
        
        % Roll-off
        cumulative_energy = cumsum(sel_power);
        roll_off_point = find(cumulative_energy >= 0.85 * cumulative_energy(end), 1, 'first');
        if ~isempty(roll_off_point)
            fast_segment_features(k, 3) = sel_freq(roll_off_point);
        else
            fast_segment_features(k, 3) = NaN;
        end
    end
end

% Combine all features and create labels
all_features = [slow_segment_features; fast_segment_features];
labels = [zeros(size(slow_segments, 1), 1); ones(size(fast_segments, 1), 1)];

% Remove NaN values if any
valid_idx = ~any(isnan(all_features), 2);
all_features = all_features(valid_idx, :);
labels = labels(valid_idx);

% Visualize the features
figure('Position', [100, 100, 1200, 400]);

% 2D plot (first two features)
subplot(1, 2, 1);
gscatter(all_features(:, 1), all_features(:, 2), labels, 'br', 'o+');
xlabel(feature_names_selected{1});
ylabel(feature_names_selected{2});
title('Feature Space Visualization (2D)');
legend('Slow', 'Fast');
grid on;

% 3D plot (all three features)
subplot(1, 2, 2);
scatter3(all_features(labels==0, 1), all_features(labels==0, 2), all_features(labels==0, 3), 'b', 'filled');
hold on;
scatter3(all_features(labels==1, 1), all_features(labels==1, 2), all_features(labels==1, 3), 'r', '+');
xlabel(feature_names_selected{1});
ylabel(feature_names_selected{2});
zlabel(feature_names_selected{3});
title('Feature Space Visualization (3D)');
legend('Slow', 'Fast');
grid on;
view(45, 30);  % Adjust viewing angle

fprintf('Classification visualization complete.\n');
fprintf('The 2D and 3D plots show how well the selected features separate slow and fast segments.\n');

% If Statistics and Machine Learning Toolbox is available, try simple classification
try
    % Create a simple classifier
    mdl = fitcdiscr(all_features, labels);
    
    % Predict and evaluate
    pred_labels = predict(mdl, all_features);
    accuracy = sum(pred_labels == labels) / length(labels);
    
    fprintf('\nLinear Discriminant Analysis Classifier:\n');
    fprintf('Classification accuracy: %.2f%%\n', accuracy * 100);
    
    % Confusion matrix
    cm = confusionmat(labels, pred_labels);
    fprintf('Confusion Matrix:\n');
    fprintf('            Predicted Slow    Predicted Fast\n');
    fprintf('Actual Slow      %d                %d\n', cm(1,1), cm(1,2));
    fprintf('Actual Fast      %d                %d\n', cm(2,1), cm(2,2));
    
    % Calculate precision, recall, etc.
    true_pos = cm(2,2);
    false_pos = cm(1,2);
    true_neg = cm(1,1);
    false_neg = cm(2,1);
    
    precision = true_pos / (true_pos + false_pos);
    recall = true_pos / (true_pos + false_neg);
    f1_score = 2 * (precision * recall) / (precision + recall);
    
    fprintf('\nClassification Metrics for "Fast" class:\n');
    fprintf('Precision: %.4f\n', precision);
    fprintf('Recall: %.4f\n', recall);
    fprintf('F1 Score: %.4f\n', f1_score);
catch
    fprintf('\nStatistics and Machine Learning Toolbox not available for classification.\n');
    fprintf('You can still analyze the feature space visualizations to see class separation.\n');
end

fprintf('\nFeature extraction and analysis in %.0f-%.0f Hz range complete.\n', ...
    freq_range_min, freq_range_max);