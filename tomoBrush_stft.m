%% Script to load 2 WAV files, find peaks, extract segments, and perform STFT Analysis
% Clear workspace and close all figures
clear all;
close all;

% Parameters for peak detection and segment extraction
target_peaks = 5; % Number of peaks to detect
smooth_window_size = 0.05; % Window size for smoothing in seconds
min_peak_distance_slow = 0.5; % Minimum distance between peaks for slow file
min_peak_distance_fast = 0.2; % Minimum distance between peaks for fast file
min_peak_prominence = 0.1; % Minimum peak prominence as fraction of max
segment_duration = 0.2; % Duration of segment to extract around each peak (in seconds)

% File names
file1 = 'TopMidBot_Slow5.wav';
file2 = 'TopMidBot_Fast5.wav';

% Process both files
processFile(file1, target_peaks, smooth_window_size, min_peak_distance_slow, min_peak_prominence, segment_duration, 'Slow');
processFile(file2, target_peaks, smooth_window_size, min_peak_distance_fast, min_peak_prominence, segment_duration, 'Fast');

function processFile(filename, target_peaks, smooth_window_size, min_peak_distance, min_peak_prominence, segment_duration, label)
    % Load audio file
    [y, Fs] = audioread(filename);
    
    % Convert to mono if stereo
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Calculate time vector
    t = (0:length(y)-1)/Fs;
    
    % Apply envelope detection
    envelope = abs(hilbert(y));
    
    % Smooth the envelope
    smooth_window = round(smooth_window_size * Fs);
    if mod(smooth_window, 2) == 0
        smooth_window = smooth_window + 1; % Make it odd for smoothing
    end
    smoothed_envelope = smoothdata(envelope, 'gaussian', smooth_window);
    
    % Find peaks
    min_peak_distance_samples = round(min_peak_distance * Fs);
    min_prominence = min_peak_prominence * max(smoothed_envelope);
    [pks, locs] = findpeaks(smoothed_envelope, 'MinPeakDistance', min_peak_distance_samples, ...
                           'MinPeakProminence', min_prominence, 'SortStr', 'descend');
    
    % Select top N peaks
    if length(pks) > target_peaks
        [~, idx] = sort(pks, 'descend');
        selected_idx = idx(1:target_peaks);
        pks = pks(selected_idx);
        locs = locs(selected_idx);
        
        % Sort by location (time) for better visualization
        [locs, sort_idx] = sort(locs);
        pks = pks(sort_idx);
    else
        fprintf('Warning: Only %d peaks found in %s (requested %d)\n', length(pks), filename, target_peaks);
    end
    
    % Display audio and detected peaks
    figure('Name', ['Audio Analysis: ' label]);
    
    % Plot waveform and peaks
    subplot(3, 1, 1);
    plot(t, y);
    hold on;
    plot(t(locs), y(locs), 'ro', 'MarkerSize', 10);
    title(['Audio Waveform with Detected Peaks: ' label]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % Plot smoothed envelope and peaks
    subplot(3, 1, 2);
    plot(t, smoothed_envelope);
    hold on;
    plot(t(locs), smoothed_envelope(locs), 'ro', 'MarkerSize', 10);
    title(['Smoothed Envelope with Detected Peaks: ' label]);
    xlabel('Time (s)');
    ylabel('Envelope Amplitude');
    grid on;
    
    % Extract segments and perform STFT analysis
    segment_samples = round(segment_duration * Fs);
    
    % Analyze each peak
    for i = 1:length(locs)
        % Calculate segment boundaries
        start_idx = max(1, locs(i) - segment_samples/2);
        end_idx = min(length(y), locs(i) + segment_samples/2);
        
        % Extract segment
        segment = y(start_idx:end_idx);
        segment_time = t(start_idx:end_idx);
        
        % Perform STFT
        figure('Name', sprintf('%s - Peak %d STFT Analysis at t=%.2fs', label, i, t(locs(i))));
        
        % Compute and display spectrogram
        window_size = min(256, length(segment)/4);
        overlap = round(window_size * 0.75);
        nfft = max(256, 2^nextpow2(window_size));
        
        % Plot segment waveform
        subplot(2, 1, 1);
        plot(segment_time, segment);
        title(sprintf('Segment Around Peak %d at t=%.2fs', i, t(locs(i))));
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
        
        % STFT/Spectrogram
        subplot(2, 1, 2);
        spectrogram(segment, window_size, overlap, nfft, Fs, 'yaxis');
        title(sprintf('STFT/Spectrogram of Peak %d', i));
        colorbar;
        
        % Alternative: stft function if available in your MATLAB version
        % [s, f, t] = stft(segment, Fs, 'Window', hamming(window_size), 'OverlapLength', overlap, 'FFTLength', nfft);
        % imagesc(t, f, 20*log10(abs(s)));
        % axis xy;
        % colorbar;
        % title(sprintf('STFT of Peak %d', i));
        % xlabel('Time (s)');
        % ylabel('Frequency (Hz)');
    end
    
    % Create a summary figure with all STFT analyses
    figure('Name', [label ' - All Peaks STFT Summary']);
    
    % Select color map
    colorMap = jet;
    
    % Determine subplot layout
    subplot_rows = ceil(sqrt(length(locs)));
    subplot_cols = ceil(length(locs) / subplot_rows);
    
    for i = 1:length(locs)
        % Calculate segment boundaries
        start_idx = max(1, locs(i) - segment_samples/2);
        end_idx = min(length(y), locs(i) + segment_samples/2);
        
        % Extract segment
        segment = y(start_idx:end_idx);
        
        % STFT parameters
        window_size = min(256, length(segment)/4);
        overlap = round(window_size * 0.75);
        nfft = max(256, 2^nextpow2(window_size));
        
        % Subplot for this peak
        subplot(subplot_rows, subplot_cols, i);
        
        % Compute STFT
        spectrogram(segment, window_size, overlap, nfft, Fs, 'yaxis');
        title(sprintf('Peak %d (t=%.2fs)', i, t(locs(i))));
        
        % Adjust colormap for better visualization
        colormap(colorMap);
    end
    
    % Adjust figure layout
    sgtitle([label ' - STFT Analysis of All Peaks']);
end