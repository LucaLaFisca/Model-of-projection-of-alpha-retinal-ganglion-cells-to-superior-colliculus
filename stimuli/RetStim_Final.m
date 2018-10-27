%% RetStimV1

%% initialization

% settings
global Set
run runSettingsEphys

% Enter rf settings here!
rfguessx = 30;    %degree
rfguessy = 15;    %degree
rf = 25;    %degree

%triggers
Set.StimB=1;
Set.TargetB=2;
Set.RewardB = 3;
dasinit(22, 2)
dasbit(Set.StimB, 0);
dasbit(Set.TargetB, 0);
dasbit(Set.RewardB, 0);
dasclearword

%cogent
cgloadlib
cgopen(Set.Screenx, Set.Screeny, 32,Set.Refresh , 1)
cgpenwid(1)
% currentTime=0;

usname = getenv('username');
Set.LogLoc = fullfile('C:','Users',usname,'Dropbox','MouseOutput','Retstim');

nts = input('Enter log name\n' ,'s'); %Enter logname

if ~isempty(nts)
    if exist([Set.LogLoc nts '.mat'], 'file')
        error('Log-file already used, please pick a different name.');
    end
else
    Set.mouse = 'Unknown';
end

%% Calculate luminance

Set.maxlum = 40;
Set.black = 0;
Set.minlum = eval([gammaconversion '(Set.black,''rgb2lum'')']);
Set.greylum = (Set.minlum+Set.maxlum)./2;
Set.white = eval([gammaconversion '(Set.maxlum,''lum2rgb'')']);
Set.grey = eval([gammaconversion '(Set.greylum,''lum2rgb'')']);

%% Parameters

%Conversion in pixels
rf = rf * Set.PixPerDeg / 2; % rf size is defined as center to edge of rf in pix
rf = ceil(rf);
rfguessx = rfguessx * Set.PixPerDeg ;    %from mouse pos
rfguessy = rfguessy * Set.PixPerDeg;    %from mouse pos
%Considering mouse position
mouseposx = (Set.ScreenWidth-Set.mouseposx) * (Set.Screenx/Set.ScreenWidth);  %from cm to px
mouseposy = Set.mouseposy * (Set.Screeny/Set.ScreenHeight); %from cm to px
rfguessx = -Set.Screenx/2 + mouseposx + rfguessx;
rfguessy = -Set.Screeny/2 + mouseposy + rfguessy;

siz_max = rf*3;
siz_min = round(rf/4);
Set.centerx = rfguessx;
Set.centery = rfguessy;
Set.siz = siz_max;

%Matrix to save each events
save_stim = struct('ID',1, 'sprite',1, 'time',0, 'centerx',Set.centerx, 'centery', Set.centery, 'size',Set.siz, 'speed',0, 'action','Init');

%% alpha RGCs parameter
smallest_diam = 4.5;   %degrees (size of stimulus that involve peak rate in alpha RGCs)
smallest_diam = smallest_diam * Set.PixPerDeg;

%% Create sprites

cgmakesprite(1,Set.siz,Set.siz,Set.grey,Set.grey,Set.grey)
cgsetsprite(1)
cgarc(0,0,Set.siz,Set.siz,0,360,[Set.white, Set.white,Set.white],'S')
cgmakesprite(2,Set.siz,Set.siz,Set.grey,Set.grey,Set.grey)
cgsetsprite(2)
cgarc(0,0,Set.siz,Set.siz,0,360,[Set.minlum, Set.minlum, Set.minlum],'S')

%Start display
cgflip(Set.grey,Set.grey,Set.grey)  %Initialisation
cgflip(Set.grey,Set.grey,Set.grey)  %colour of background
cgsetsprite(0)  %Next drawings on default window
cgflip(Set.grey,Set.grey,Set.grey)
cgflip(Set.grey,Set.grey,Set.grey)
disp('Press any button to start')
pause

time = tic;

%% Create trial structure

count = 1;
stimB_count = 0;
ID_count = 0;
siz = Set.siz;
stim_dur = 2; %sec
iti=tic; %initialize iti
ititime = 1;
fixed_time = 0.5; %fix the sprite in the center before each trial

%clear keyboard
[kd,kp] = cgkeymap;

%%ID definition:
%1->3: up/down
%4->6: left/right
%7->9: increase/decrease
%1->9: white
%10->18: black
%1->18: size = center_diam
%19->36: size = medium_diam
%37->54: size = biggest_diam
num_ID = 54;

%initialize random trial display
ID_vect = randperm(num_ID);

%initialise random sizes sequence
biggest_diam = rf*3;
medium_diam = (smallest_diam+biggest_diam)/2;

ESC=0;

%Run until ESC is pressed
while ~ESC
    %     count
    targetB_count = 0;
    stimbitsent = 0;
    
    %Set the stimulus
    ID = ID_vect(1);
    siz_choice = ceil(ID/18);
    switch siz_choice
        case 1
            easy_ID1 = ID;
            siz = smallest_diam;
        case 2
            easy_ID1 = ID - 18;
            siz = medium_diam;
        case 3
            easy_ID1 = ID - 2*18;
            siz = biggest_diam;
    end
    sprite = ceil(easy_ID1/9);
    easy_ID2 = mod(easy_ID1,9);
    if ~easy_ID2
        easy_ID2 = 9;
    end
    choice = ceil(easy_ID2/3);
     speed = (1+(easy_ID2 - 3*(choice-1)))^2;
    
    %Run
    stimB_count = stimB_count + 1;
    fprintf('Trial %i \n',stimB_count);
    dasword(ID);
    centery = Set.centery;
    centerx = Set.centerx;
    if toc(iti)<ititime
        pause(ititime - toc(iti));
    end;
    cgdrawsprite(sprite,centerx,centery,siz,siz);
    cgflip(Set.grey,Set.grey,Set.grey)
    dasbit(Set.RewardB, 1);
    save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',0, 'action','Appear');
    count = count+1;
    pause(fixed_time)
    dasbit(Set.RewardB, 0);
    switch choice
        case 1  %Upwards/Downwards
            motion_dist = speed; 
            up_first = randi([0 1]);
            prev_count = count;
            disptime = tic;
            while toc(disptime) < stim_dur
                pause(0.001)
                if up_first
                    dasbit(Set.TargetB,0);
                    while and(centery <= Set.centery + rf - motion_dist, ...
                            toc(disptime) < stim_dur)    %up
                        centery = centery + speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Upwards');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count;
                    end
                    up_first = 0;
                    
                else
                    dasbit(Set.TargetB,0);
                    while and(centery >= Set.centery - rf + motion_dist, ...
                            toc(disptime) < stim_dur)   %down
                        centery = centery - speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Downwards');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count;
                    end
                    up_first = 1;
                end
            end
        case 2  %left/right
            motion_dist = speed;
            left_first = randi([0 1]);
            prev_count = count;
            disptime = tic;
            while toc(disptime) < stim_dur
                pause(0.001)
                if left_first
                    dasbit(Set.TargetB,0);
                    while and(centerx >= Set.centerx - rf + motion_dist, ...
                            toc(disptime) < stim_dur)   %left
                        centerx = centerx - speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Left');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count;
                    end
                    left_first = 0;
                else
                    dasbit(Set.TargetB,0);
                    while and(centerx <= Set.centerx + rf - motion_dist, ...
                            toc(disptime) < stim_dur)   %right
                        centerx = centerx + speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Right');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count;
                    end
                    left_first = 1;
                end
            end
        case 3  %increase/decrease
            increase_first = randi([0 1]);
            disptime = tic;
            prev_count = count;
            while toc(disptime) < stim_dur
                pause(0.001)
                if increase_first
                    dasbit(Set.TargetB,0);
                    while and(siz <= siz_max - speed,  ...
                            toc(disptime) < stim_dur)     %increase
                        siz = siz + speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Increase');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count; 
                    end
                    increase_first = 0;
                else
                    dasbit(Set.TargetB,0);
                    while and(siz >= siz_min + speed, ...
                            toc(disptime) < stim_dur)     %decrease
                        siz = siz - speed;
                        cgdrawsprite(sprite,centerx,centery,siz,siz);
                        cgflip(Set.grey,Set.grey,Set.grey)
                        if ~stimbitsent
                            dasbit(Set.StimB,1);
                            stimbitsent=1;
                        end
                        save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',speed, 'action','Decrease');
                        count = count+1;
                    end
                    if not(prev_count == count)
                        targetB_count = targetB_count + 1;
                        dasbit(Set.TargetB,1);
                        prev_count = count;
                    end
                    increase_first = 1;
                end
            end
            %         case 4  %change sprite
            %             sprite = mod(sprite,2) + 1;
            %             %gray background
            %             cgflip(Set.grey,Set.grey,Set.grey)
            %             cgflip(Set.grey,Set.grey,Set.grey)
            %             save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',0, 'action','Disappear');
            %             pause(stim_dur/2)
            %             targetB_count = targetB_count + 1;
            %             dasbit(Set.TargetB,1);
            %             %changed stimulus
            %             cgdrawsprite(sprite,centerx,centery,siz,siz);
            %             cgflip(Set.grey,Set.grey,Set.grey)
            %             if ~stimbitsent
            %                 dasbit(Set.StimB,1);
            %                 stimbitsent=1;
            %             end
            %             save_stim(count) = struct('ID',ID, 'sprite',sprite, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',0, 'action','Change');
            %             pause(stim_dur/2)
            %             count = count+1;
            %             dasbit(Set.TargetB,0);
    end
    cgflip(Set.grey,Set.grey,Set.grey)
    dasbit(Set.RewardB, 1);
    save_stim(count) = struct('ID',ID, 'sprite',0, 'time',toc(time), 'centerx',centerx, 'centery', centery, 'size',siz, 'speed',0, 'action','Disappear');
    count = count+1;
    iti = tic;
    dasclearword;
    dasbit(Set.StimB,0);
    dasbit(Set.TargetB,0);
    dasbit(Set.RewardB, 0);
    ID_count = ID_count + 1;
    if mod(ID_count,num_ID)==0;
        fprintf('all cases done %i times \n', ID_count/num_ID);
    end
    
    ID_vect(1)=[];
    
    %new random order after each stimulus has been displayed once
    if isempty(ID_vect)
        disp('Randomizing...')
        ID_vect = randperm(num_ID);
    end
    
    %check for ESC press
    [kd,kp] = cgkeymap;
    if length(find(kp)) == 1
        if find(kp) == 1;
            ESC = 1;
        end
    end
end

toc(time)
cgfreesprite([1 2])
if ~isempty(nts)
    try
        save(fullfile(Set.LogLoc, nts),'save_stim');
    catch me
        disp(me);
        save(nts,'save_stim')
    end
end
cgshut
clear all 
