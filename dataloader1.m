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

%% Main Script
file_name = 'A.wav';
dir_path = 'data2'; 
full_path = fullfile(dir_path, file_name);
[audio, fs] = audioread(full_path);
if size(audio,2) > 1
    audio = audio(:, 1); % Convert to mono if stereo
end
audio = squeeze(audio);
t = (0:length(audio)-1) / fs;

%% Segmentation + feature check
segment_duration = 0.15; % s
freq_res = fs / 4096;
freq_range = [0.2, 10000]; % kHz

% peak finder
MinPeakHeight = 0.05;
[peak_values, peak_indices] = findpeaks(audio, 'MinPeakHeight',MinPeakHeight, 'MinPeakDistance', 5000);
segments = zeros(length(peak_indices), segment_duration*fs*2+1);
n_seg = length(peak_indices);

figure
for i = 1:n_seg
    subplot(n_seg, 3, i*3-2)
    seg = audio(peak_indices(i)-segment_duration*fs:peak_indices(i)+segment_duration*fs);
    segments(i,:) = seg;
    plot(seg)
    subplot(n_seg, 3, i*3-1)
    pspectrum(seg, fs, 'FrequencyResolution', freq_res)
    xlim(freq_range/1000)
    subplot(n_seg, 3, i*3)
    plot_powerSpectrum(seg, fs, 1);
end

%% create dataset
% init
[Pxx, f] = pspectrum(segments(1,:), fs, 'FrequencyResolution', freq_res);
freq_mask = f >= freq_range(1) & f <= freq_range(2);
X = zeros(n_seg, sum(freq_mask));

% retrieve all
for i = 1:n_seg
    [Pxx, f] = pspectrum(segments(i,:), fs, 'FrequencyResolution', freq_res);
    X(i,:) = Pxx(freq_mask);
end

X = log10(X + eps);

% get label - Modified for braille pattern
[~, name_only, ~] = fileparts(file_name);
letter = name_only; % Just the letter (A, B, C, etc.)
braille_pattern = get_braille_pattern(letter);
class_label = braille_to_class(braille_pattern);
y = repmat(class_label, n_seg, 1);

% check
fprintf("Feature matrix size: %d samples × %d freq bins\n", size(X,1), size(X,2));
fprintf("Letter: %s, Braille pattern: [%s], Class: %d\n", letter, ...
    num2str(braille_pattern), class_label);

%% create all dataset
dir_path = 'data2';  % Changed to data2
file_list = dir(fullfile(dir_path, '*.wav'));
X_all = [];
y_all = [];

for file = file_list'
    file_name = file.name;
    full_path = fullfile(dir_path, file_name);
    [audio, fs] = audioread(full_path);
    if size(audio, 2) > 1
        audio = audio(:,1); % convert to mono
    end
    
    [peak_values, peak_indices] = findpeaks(audio, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 5000);
    n_seg = length(peak_indices);
    
    for i = 1:n_seg
        idx_start = peak_indices(i) - segment_duration*fs;
        idx_end = peak_indices(i) + segment_duration*fs;
        if idx_start < 1 || idx_end > length(audio)
            continue;
        end
        seg = audio(idx_start:idx_end);
        [Pxx, f] = pspectrum(seg, fs, 'FrequencyResolution', fs/512);
        freq_mask = f >= freq_range(1) & f <= freq_range(2);
        spectrum_feature = Pxx(freq_mask)';
        X_all = [X_all; log10(spectrum_feature + eps)];
    end
    
    % Modified labeling for braille patterns
    [~, name_only, ~] = fileparts(file_name);
    letter = name_only;
    braille_pattern = get_braille_pattern(letter);
    class_label = braille_to_class(braille_pattern);
    
    % Add labels for all segments from this file
    y_labels = repmat(class_label, size(X_all,1) - length(y_all), 1);
    y_all = [y_all; y_labels];
end

% Display class distribution
unique_classes = unique(y_all);
fprintf('\nClass distribution:\n');
for i = 1:length(unique_classes)
    count = sum(y_all == unique_classes(i));
    fprintf('Class %d: %d samples\n', unique_classes(i), count);
end

fprintf('\nTotal unique classes: %d\n', length(unique_classes));

save('ML_output/Xy_data1.mat', 'X_all', 'y_all');

%% Optional: Create a lookup table for class to letter mapping
class_to_letter = containers.Map('KeyType', 'int32', 'ValueType', 'char');
letter_to_class = containers.Map();

for letter = 'A':'Z'
    braille_pattern = get_braille_pattern(letter);
    class_label = braille_to_class(braille_pattern);
    class_to_letter(class_label) = letter;
    letter_to_class(letter) = class_label;
end

save('ML_output/braille_mapping.mat', 'class_to_letter', 'letter_to_class');