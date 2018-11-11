%% Load EVENT for RFmap
Word = EVENT.Trials.word;
[row, ~] = find(isnan(Word));
wrong = sort(unique(row));
% %% Load Mat, cleanMua and px for RFmap
% RFmapMua = cleanMua;
% RFmapMat = Mat;
% RFmapPx = px*1000; %Convert to ms
% 
% %% 
% clear cleanMua
% clear Mat
% clear px
% %% Load Mat, cleanMua and px for RetStim
% RetstimMua1 = cleanMua;
% % RetstimMua2 = cleanMua2;
% RetstimMat = struct2table(save_stim);
% RetstimPx = px*1000; %Convert to ms
% 
% for i = 1:32
%     RFmapMua{i,1}(:,wrong) = [];
%     nanDetect1 = find(isnan(RetstimMua1{i,1}));
% %     nanDetect2 = find(isnan(RetstimMua2{i,1}));
% end
% 
% % if(and(isempty(nanDetect1),isempty(nanDetect2)))
% if isempty(nanDetect1)
%     disp('No NaN detected')
% else
%     disp('WARNING !!! NaN detected !')
%     pause()
% end
% 
% if ~(size(RFmapMat,1) == size(RFmapMua{1,1},2))
%     disp('Mismatch between data and log file!')
%     disp('Solve it manually thanks to "wrong" vector and "Word" array')
%     pause()
% end

%% RFmap
dt = .5; %ms
RFmapMua = cleanMua;
RFmapMat = Mat;
RFmapPx = px*1000; %Convert to ms
pointsPerTrial = (RFmapPx(end)-RFmapPx(1))/dt;
RFmapStep = floor(size(RFmapMua{1,1},1)/pointsPerTrial);
RFmapPx = RFmapPx(1:RFmapStep:end);
% RFmapMUA = zeros(size(RFmapMua{1,1},1),size(RFmapMua{1,1},2),...
%     size(RFmapMua,1));


for i = 1:size(RFmapMua,1)
    RFmapMUA(:,:,i) = RFmapMua{i,1}(1:RFmapStep:end,:);
%     RetstimMUA(:,:,i) = horzcat(RetstimMua1{i,1}(1:RetstimStep:end,:),RetstimMua2{i,1}(1:RetstimStep:end,:));
end

% RFmapMUA(1,:,:) = [];   %Because the first value of each trial is weird due to filtering
% RFmapPx(1) = [];
%%
save('Recordings\Session-5\RFmapPx.mat','RFmapPx')
save('Recordings\Session-5\RFmapMat.mat','RFmapMat')
save('Recordings\Session-5\RFmapMUA.mat','RFmapMUA')


%% RetStim
dt=.5;
RetstimMua = cleanMua;
RetstimMat = struct2table(save_stim);
RetstimPx = px*1000; %Convert to ms
pointsPerTrial = (RetstimPx(end)-RetstimPx(1))/dt;
RetstimStep = floor(size(RetstimMua{1,1},1)/pointsPerTrial);
RetstimPx = RetstimPx(1:RetstimStep:end);

% RetstimMUA = zeros(size(RetstimMua1{1,1},1),...
%     size(RetstimMua1{1,1},2)+size(RetstimMua1{1,1},2),size(RetstimMua1,1));

for i = 1:size(RetstimMua,1)
    RetstimMUA(:,:,i) = RetstimMua{i,1}(1:RetstimStep:end,:);
%     RetstimMUA(:,:,i) = horzcat(RetstimMua1{i,1}(1:RetstimStep:end,:),RetstimMua2{i,1}(1:RetstimStep:end,:));
end


%%

save('Recordings\Session-5\RetStimPx.mat','RetstimPx')
save('Recordings\Session-5\RetStimMat.mat','RetstimMat')
save('Recordings\Session-5\RetStimMUA.mat','RetstimMUA')