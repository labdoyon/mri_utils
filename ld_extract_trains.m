% Detection of grouped and isolated spindles
% Extract variables such as ISI, number of trains, etc.
% Write markers (grouped, isolated) into vmrk

clear all 
clc

mainDir = 'D:\Recherche\CRIUGM\Sleep_Reconsolidation\Analysis_Night1\Data\spindles\vmrk_grouped_isolated\';

spindleFolder = fullfile(mainDir, 'spindles');
outputFolder = fullfile(mainDir, 'output');

% Electrode of interest
o_vars.currentElectrode = 'Pz';

% SleepStage 
o_vars.sleepStage = '2'; %('2' for NREM2, '3' for NREM3 and '23' for NREM2 and NREM3)
if strcmp(o_vars.sleepStage,'2')
    indexSleepStage = '3'; % Hardcoded to check
    disp('Sleep stage: NREM2') 
elseif strcmp(o_vars.sleepStage,'3')
    indexSleepStage = '4'; % Hardcoded to check
    disp('Sleep stage: NREM3') 
elseif strcmp(o_vars.sleepStage,'23')
    indexSleepStage = '3:4'; % Hardcoded to check
	disp('Sleep stage: NREM2 and NREM3') 
end

% Subjects
allFiles = dir(fullfile(spindleFolder, '*.mat'));

% Empty structures
o_vars.subjects = cell(1,length(allFiles));

% Save subject Names
o_vars.subjectNames = cell(1,length(allFiles));

% Detection trains parameters
o_vars.ISI = 6;
o_vars.nbMinSpPerTrain = 2;

% Loop over all subjects / files
for nFile=1:length(allFiles)
    o_vars.subjects{nFile}.filename = allFiles(nFile).name;
    
    o_vars.subjects{nFile}.subjectName = strsplit(allFiles(nFile).name(3:end),'_Day');
    o_vars.subjects{nFile}.subjectName = o_vars.subjects{nFile}.subjectName{1};
    
    o_vars.subjectNames{nFile} = o_vars.subjects{nFile}.subjectName;
    
    % Display file name
    disp('%%%%%%%%%%')
    disp(['Filename: ', allFiles(nFile).name])
    disp(['Subject: ', o_vars.subjects{nFile}.subjectName])
    
    % LOAD Spindles
    %disp('Load spindles')
    i_SS = load(fullfile(spindleFolder, allFiles(nFile).name), 'SS');
    i_SS = i_SS.SS;
    
    % LOAD Info
    i_Info = load(fullfile(spindleFolder, allFiles(nFile).name), 'Info');
    i_Info = i_Info.Info;
    
    indexElectrode = find(strcmp({i_Info.Electrodes.labels}, o_vars.currentElectrode)==1);
    
    % Show all spindles
    disp(['Found ', num2str(length(i_SS)) ' spindles'])
    
    % Found indexes with correct electrodes and correct scoring
    scoring = reshape([i_SS.scoring], length(i_Info.Electrodes), length(i_SS))';
    refStart = reshape([i_SS.Ref_Start], length(i_Info.Electrodes), length(i_SS))';
    
    % Found indexes
    indexSpFound = find(scoring(:,indexElectrode) == str2double(o_vars.sleepStage));
    
    % Get spindles
    currSS = i_SS(indexSpFound);
    
    % Show spindles on specific electrode and scoring
    disp(['Found ' num2str(length(currSS)), ' spindles on ' o_vars.currentElectrode])
        
    % Get all beginings
    currSpStarts = sort(refStart(indexSpFound, indexElectrode));

    %%
    % Compute groups
    spGroup = {};
    spNoGroup = {};
    
    lastStart = currSpStarts(1); % first spindle of the first group
    groupIndexes = 1; % first group
    
    for nSp=2:length(currSpStarts)
        if (currSpStarts(nSp)-lastStart)/i_Info.Recording.sRate < o_vars.ISI % New group
            groupIndexes(end+1) = nSp;
            lastStart = currSpStarts(nSp);
        elseif length(groupIndexes) < o_vars.nbMinSpPerTrain % Valid spindle for group
            spNoGroup{end+1} = groupIndexes;
            lastStart = currSpStarts(nSp);
            groupIndexes = nSp;
        else % End of the group
            spGroup{end+1} = groupIndexes;
            lastStart = currSpStarts(nSp);
            groupIndexes = nSp;
        end
    end
    
    if length(groupIndexes) < o_vars.nbMinSpPerTrain % Last group
        spNoGroup{end+1} = groupIndexes;
    else
        spGroup{end+1} = groupIndexes;
    end
    
    % Init interTrain
    o_vars.subjects{nFile}.interTrain = 0;
    
    %%
    % Store grouped and isolated spindles
    for nSpG=1:length(spGroup)
        currGroup = spGroup{nSpG};
        startsGroup = [];
        for currSp=1:length(currGroup)
            startsGroup(end+1) =  currSS(currGroup(currSp)).Ref_Start(indexElectrode);
        end
        
        startsGroup = sort(startsGroup);
        
        diffStartsGroup = diff(startsGroup)/i_Info.Recording.sRate;
        
        try
            o_vars.subjects{nFile}.spGroup(end+1).min = min(diffStartsGroup);
        catch
            o_vars.subjects{nFile}.spGroup(1).min = 0;
            o_vars.subjects{nFile}.spGroup(end).min = min(diffStartsGroup);
        end
        
        o_vars.subjects{nFile}.spGroup(end).max = max(diffStartsGroup);
        o_vars.subjects{nFile}.spGroup(end).mean = mean(diffStartsGroup);
        o_vars.subjects{nFile}.spGroup(end).median = median(diffStartsGroup);
        o_vars.subjects{nFile}.spGroup(end).length = sum(diffStartsGroup)+currSS(currGroup(currSp)).Ref_Length(indexElectrode)/i_Info.Recording.sRate;
        o_vars.subjects{nFile}.spGroup(end).nb = length(currGroup);
        
        if nSpG > 1 
            o_vars.subjects{nFile}.interTrain(end+1) = (startsGroup(1)-endGroup)/i_Info.Recording.sRate;
        end    
        
        endGroup = startsGroup(end) + currSS(currGroup(currSp)).Ref_Length(indexElectrode);            
    end
    
    % Remove first dummy interTrain
    o_vars.subjects{nFile}.interTrain(1) = [];
    
    o_vars.subjects{nFile}.spindlesGrouped.indexes = cell2mat(spGroup);
    o_vars.subjects{nFile}.spindlesGrouped.nb = sum([o_vars.subjects{nFile}.spGroup(:).nb]);
    o_vars.subjects{nFile}.spindlesNotGrouped.indexes = cell2mat(spNoGroup);
    o_vars.subjects{nFile}.spindlesNotGrouped.nb = length(spNoGroup);
    
    disp(['Number of trains : ' num2str(length(spGroup)) ' with a total of '  num2str(sum([o_vars.subjects{nFile}.spGroup(:).nb])) ' spindles'])
    disp(['Not grouped Spindles : ' num2str(length(spNoGroup))])
    
    %% vmrk loop (start)
    % Save to VMRK Uncomment
    % Structure of a marker
    
    o_vars.subjects{nFile}.o_vmrkFileName = strrep(o_vars.subjects{nFile}.filename,'.mat','.vmrk'); 
    
    Marker = struct('type',{},'description',{},'position',{},'length',{},'channel',{});
    
    % Isolated spindles 
    for nSp=1:length(o_vars.subjects{nFile}.spindlesNotGrouped.indexes)
        currIndex = o_vars.subjects{nFile}.spindlesNotGrouped.indexes(nSp);
        newmark = struct('type','SpNotGrouped', ...
                        'description',['Sp_' currSS(nSp).Ref_TypeName{indexElectrode} '_' o_vars.currentElectrode '_NREM' num2str(currSS(nSp).scoring(indexElectrode))], ...
                        'position',currSS(currIndex).Ref_Start(indexElectrode), ...
                        'length',currSS(currIndex).Ref_Length(indexElectrode), ...
                        'channel',indexElectrode); 
        Marker = [Marker newmark];                    
    end
    
    % Grouped spindles
    for nSp=1:length(o_vars.subjects{nFile}.spindlesGrouped.indexes)
        currIndex = o_vars.subjects{nFile}.spindlesGrouped.indexes(nSp);
        newmark = struct('type','SpGrouped', ...
                        'description',['Sp_' currSS(nSp).Ref_TypeName{indexElectrode} '_' o_vars.currentElectrode '_NREM' num2str(currSS(nSp).scoring(indexElectrode))], ...
                        'position',currSS(currIndex).Ref_Start(indexElectrode), ...
                        'length',currSS(currIndex).Ref_Length(indexElectrode), ...
                        'channel',indexElectrode); 
        Marker = [Marker newmark];                    
    end
    
    fid = fopen(o_vars.subjects{nFile}.o_vmrkFileName ,'w');
    fprintf(fid,'%s\n\n','Brain Vision Data Exchange Marker File, Version 2.0'...
        ,'; Data created from history path:'...
        ,'; The channel numbers are related to the channels in the exported file.');

    fprintf(fid,'%s\n','[Common Infos]', ...
               'Codepage=UTF-8');
    
    datFileName = strsplit(o_vars.subjects{nFile}.o_vmrkFileName,filesep);
    datFileName = datFileName{end};
    fprintf(fid,'%s\n\n',['DataFile=',datFileName]);

    fprintf(fid,'%s\n','[Marker Infos]'...
    ,'; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,'...
    ,'; <Size in data points>, <Channel number (0 = marker is related to all channels)>'...
    ,'; Fields are delimited by commas, some fields might be omitted (empty).'...
    ,'; Commas in type or description text are coded as');

    for i = 1:length(Marker)
        fprintf(fid,'%s\n',['Mk' num2str(i) '=' Marker(i).type ',' Marker(i).description ',' num2str(Marker(i).position) ',' num2str(Marker(i).length) ',' num2str(Marker(i).channel)]);
    end

    %%% vmrk loop (end)
end

%% Merge o_vars stats 

[D,~,X] = unique(o_vars.subjectNames(:));
Y = hist(X,unique(X));
Z = struct('name',D,'freq',num2cell(Y(:)));

spindlesAnalysis = o_vars;
spindlesAnalysis.subjects = cell(1,length(D));

spindlesAnalysis.subjectNames = D;

for nSub=1:length(D)
    indexes = find(strcmp(o_vars.subjectNames, D(nSub)));

    spindlesAnalysis.subjects{nSub} = o_vars.subjects{indexes(1)};
    spindlesAnalysis.subjects{nSub}.spindlesGrouped = spindlesAnalysis.subjects{nSub}.spindlesGrouped.nb;
    spindlesAnalysis.subjects{nSub}.spindlesNotGrouped = spindlesAnalysis.subjects{nSub}.spindlesNotGrouped.nb;

    spindlesAnalysis.subjects{nSub}.filename = {o_vars.subjects{nSub}.filename};
    
    if Z(nSub).freq > 1      
        for nIndex = 2:length(indexes)
            
            spindlesAnalysis.subjects{nSub}.filename{end+1} = o_vars.subjects{indexes(nIndex)}.filename;
            
            % Basic info spindles
            spindlesAnalysis.subjects{nSub}.spindlesGrouped = spindlesAnalysis.subjects{nSub}.spindlesGrouped + o_vars.subjects{indexes(nIndex)}.spindlesGrouped.nb;
            spindlesAnalysis.subjects{nSub}.spindlesNotGrouped = spindlesAnalysis.subjects{nSub}.spindlesNotGrouped + o_vars.subjects{indexes(nIndex)}.spindlesNotGrouped.nb;
            spindlesAnalysis.subjects{nSub}.interTrain = [spindlesAnalysis.subjects{nSub}.interTrain,  o_vars.subjects{indexes(nIndex)}.interTrain];
            
            % Infos spindles
            spindlesAnalysis.subjects{nSub}.spGroup = [spindlesAnalysis.subjects{nSub}.spGroup, o_vars.subjects{indexes(nIndex)}.spGroup];
        end
    end

    save('out_analysis.mat', 'spindlesAnalysis');    
end