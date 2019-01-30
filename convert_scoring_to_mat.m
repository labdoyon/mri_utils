% Extract scoring to mat file - Benoit Adam
% 
% Arnaud Boré - arnaud.bore@gmail.com

%% Try to replicate markers structure from old scoring Laura Ray

% length epoch in seconds
% D.other.CRC.score{1,1} =

% Scoring
% markers.D.other.CRC.score{1,1} = 

scoringFolder = '/home/bore/p/deconap/eeg_analysis';

%% READ Files
allFiles = dir(fullfile(scoringFolder, '*.csv'));

for nFile=1:length(allFiles)
    
    iPath = fullfile(allFiles(nFile).folder, allFiles(nFile).name);
    
    splitFilename = strsplit(allFiles(nFile).name,'_');
    subject = splitFilename{2};
    subjectSess = splitFilename{4};
    
    scoringTable = readtable(iPath);

    D.other.CRC.score{1,1} = scoringTable.Stade';
    D.other.CRC.score{2,1} = {};

    epoch = strsplit(char(scoringTable.Dur_e_Sec_(1)),',');
    epoch = str2double(epoch(1));

    D.other.CRC.score{3,1} = epoch;

    save(['DeCoNap_' subject '_Nap_' subjectSess '_scoring.mat'],'D');
end
