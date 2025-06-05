function braille_pattern = get_braille_pattern(letter)
    % Braille dot positions:
    % 1 4
    % 2 5
    % 3 6
    
    braille_map = containers.Map();
    braille_map('A') = [1];
    braille_map('B') = [1,2];
    braille_map('C') = [1,4];
    braille_map('D') = [1,4,5];
    braille_map('E') = [1,5];
    braille_map('F') = [1,2,4];
    braille_map('G') = [1,2,4,5];
    braille_map('H') = [1,2,5];
    braille_map('I') = [2,4];
    braille_map('J') = [2,4,5];
    braille_map('K') = [1,3];
    braille_map('L') = [1,2,3];
    braille_map('M') = [1,3,4];
    braille_map('N') = [1,3,4,5];
    braille_map('O') = [1,3,5];
    braille_map('P') = [1,2,3,4];
    braille_map('Q') = [1,2,3,4,5];
    braille_map('R') = [1,2,3,5];
    braille_map('S') = [2,3,4];
    braille_map('T') = [2,3,4,5];
    braille_map('U') = [1,3,6];
    braille_map('V') = [1,2,3,6];
    braille_map('W') = [2,4,5,6];
    braille_map('X') = [1,3,4,6];
    braille_map('Y') = [1,3,4,5,6];
    braille_map('Z') = [1,3,5,6];
    
    if isKey(braille_map, upper(letter))
        braille_pattern = braille_map(upper(letter));
    else
        braille_pattern = [];
    end
end

function class_label = braille_to_class(braille_pattern)
    % Each pattern gets a unique number based on binary representation
    class_vector = zeros(1, 6);
    class_vector(braille_pattern) = 1;
    % Convert binary pattern to decimal -> unique class id
    class_label = sum(class_vector .* (2.^(0:5)));
end
%% === 1) Extract + label every segment, storing time‐domain “segments_all” ===
% (This replaces your “for file = file_list” loop, with an added container
%  segments_all that keeps the raw segment waveforms. After this loop,
%  segments_all(i,:) is the i-th segment’s waveform, and y_all(i) is its label.)

dir_path = 'data2';  
file_list = dir(fullfile(dir_path, '*.wav'));
MinPeakHeight = 0.05;
segment_duration = 0.15; % seconds

segments_all = {};  % we’ll store each segment’s raw waveform here
y_all = [];         % matching class labels

for file = file_list'
    % Read in audio
    [audio, fs] = audioread(fullfile(dir_path, file.name));
    if size(audio,2) > 1
        audio = audio(:,1);
    end
    
    % Find peaks and segment around each peak
    [~, peak_indices] = findpeaks(audio, 'MinPeakHeight', MinPeakHeight, ...
                                  'MinPeakDistance', 5000);
    n_seg = length(peak_indices);

    for i = 1:n_seg
        idx0 = peak_indices(i);
        idx_start = idx0 - round(segment_duration*fs);
        idx_end   = idx0 + round(segment_duration*fs);
        if idx_start < 1 || idx_end > numel(audio)
            continue; 
        end
        
        seg = audio(idx_start:idx_end);
        segments_all{end+1,1} = seg(:)';  % store it as a row‐vector
    end
    
    % Compute this file’s class label (based on filename/letter)
    [~, base, ~] = fileparts(file.name);
    braille_pattern = get_braille_pattern(base);
    cls = braille_to_class(braille_pattern);
    
    % Assign that same class to all segments just extracted
    labels = repmat(cls, sum(~cellfun(@isempty, segments_all((end-n_seg+1):end))), 1);
    y_all = [y_all; labels];
end

% Now we have:
%   segments_all = N×1 cell array, each row is a segment waveform
%   y_all         = N×1 vector of integer class IDs

%% === 2) Loop over each unique class and plot its segments ===

unique_classes = unique(y_all);

for k = 1:numel(unique_classes)
    cls = unique_classes(k);
    idx_cls = find(y_all == cls);
    n_cls  = numel(idx_cls);
    
    % Create a new figure for this class
    figure('Name', sprintf('Class %d: %d segments', cls, n_cls), ...
           'NumberTitle','off');
    
    % Decide grid size for subplots (e.g. up to 5 columns)
    ncols = min(5, n_cls);
    nrows = ceil(n_cls / ncols);
    
    for j = 1:n_cls
        seg_wave = segments_all{idx_cls(j)};
        t = (0:length(seg_wave)-1) / fs;  % time axis (seconds)
        
        subplot(nrows, ncols, j);
        plot(t, seg_wave, 'LineWidth', 1);
        xlabel('Time (s)'); 
        ylabel('Amplitude');
        title(sprintf('Seg %d', j));
        axis tight;
    end
    
    % If you also want the frequency‐domain next to each time plot, 
    % you could replace the “plot(t,seg_wave)” line with something like:
    %
    %   subplot(nrows, 2*ncols, 2*j-1);
    %   plot(t, seg_wave); xlabel('Time (s)'); title(sprintf('Time %d',j));
    %
    %   subplot(nrows, 2*ncols, 2*j);
    %   [Pxx,f] = pspectrum(seg_wave, fs, 'FrequencyResolution', freq_res);
    %   plot(f, 10*log10(Pxx)); xlim(freq_range); 
    %   xlabel('Freq (Hz)'); title(sprintf('PSD %d',j));
    %
    % and adjust the grid dimensions accordingly (2*ncols across).
    
    % Pause after each figure so you can inspect before moving on:
    fprintf('Plotted %d segments for class %d. Press any key to continue...\n', n_cls, cls);
    pause;
end