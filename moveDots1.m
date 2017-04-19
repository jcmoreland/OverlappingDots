% moveDots1.m
% Display showing two fields of overlapping dots. One on each side of the
% display.

clear all
meta.expdir = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\FeatureBasedAttention\OverlappingDots';

meta.subj = 'KM';
meta.expName = 'moveDots1';
meta.date = datestr(now, 'dd-mmm-yyyy');            % date
meta.scriptDate = '04/18/17';
meta.computer = 0;  % 0 = laptop, 1 = lab
meta.debugging = 1;
meta.initialSeed = ClockRandSeed;

fprintf('\nSubject: %s\n', meta.subj);
fprintf('\nExperiment: %s\n', meta.expName);

if ~meta.computer
    p.computerName = 'DellXPS13';
else
    p.computerName = 'coombs';
end
%% Folder Paths

meta.dataDr = fullfile(meta.expdir,meta.subj, 'Data');
if ~isdir(meta.dataDr)
    mkdir(meta.dataDr)
end
meta.noiseFolder = fullfile(meta.expdir,'movies/');

% Set up save files
session = '_session_';

meta.saveFileName = [meta.expName,'_',meta.subj,session];
meta.sessionFile = fullfile(meta.dataDr,meta.saveFileName);

meta.blockTmpFileName = [meta.expName, '_', meta.subj, '_tmp_'];
meta.BlockTmpFile = fullfile(meta.dataDr,meta.blockTmpFileName); % This format will be used to add blocks onto the end of.

%% Condition structure
if ~meta.debugging
    cond.nTrials = 30;                  % Number of trials per block
    cond.blockatt = [1,2,1,2,1,2,1,2];          % Attention conditions 1==same 2==diferent
    cond.blockleftdir = [1,1,2,2,1,1,2,2];      % Motion direction of attended left (if same then right is same direction etc)
    cond.nBlocks = length(cond.blockatt);% Number of blocks
else
    cond.nTrials = 2;                   % Number of trials per block
    cond.blockatt = [2];          % Attention conditions 1==same 2==diferent
    cond.blockleftdir = [1];      % Motion direction of attended left (if same then right is same direction etc)
    cond.nBlocks = length(cond.blockatt);% Number of blocks
end

%% Response Keys

buttons = getKeyAssignment(p.keyboard);
p.buttons = buttons;

%% Display Parameters:
display.screenNum = max(Screen('Screens'));
display.dist = 50;      % cm office 50
display.width = 37;     % cm office 42?
% display.bkColor = [128 128 128];
display.bkColor = [255 255 255];
display.textColor = round(.75*[255 255 255]);
display.fixColor = [0 0 0];
display.fixOffset = [0 0];
display.fontSize = 30;
display.fontName = 'Courier';

Screen('Preference', 'SkipSyncTests',display.screenNum);
Screen('Preference','TextRenderer',1);

waitframes = 1; % flip every other frame

%% Open Screen

[display.windowPtr, rect] = Screen('OpenWindow',display.screenNum,display.bkColor);
topPriorityLevel = MaxPriority(display.windowPtr);
display.ifi = Screen('GetFlipInterval',display.windowPtr);
display.adj = (waitframes - 0.5)*display.ifi;
display.frameRate = 1/display.ifi; %Hz
display.resolution = rect([3,4]);
display.center = display.resolution/2;
Screen('TextSize', display.windowPtr, display.fontSize);
Screen('TextFont', display.windowPtr, display.fontName);
Screen('TextStyle', display.windowPtr, 0);

scr.monParamFileName = sprintf('expMonitorParams_%s.mat',p.computerName);
if ~exist(scr.monParamFileName)
    scr.monParamFileName = 'expMonitorParams_default';
    fprintf(1,'\n\n\n(prepScreen) !!!!!! Warning !!!!! \n(prepScreen)\t No monitor params file for this computer! Using defualt!!!\n\n');
end
load(scr.monParamFileName);
if meta.computer
    display.cl = CCGetCalibration(monParams.calibFile);
    gammaInverse = CCMakeInverseGamma(display.cl);
    BackupCluts; % This command restores the CLUTs with sca later.
    Screen('LoadNormalizedGammaTable', display.windowPtr, gammaInverse);
    fprintf('\n(prepScreen) Loaded calibration file normalized gamma table: %s\n',monParams.calibFile);
    fprintf(1,'--------------------------------------------------------------\n');
end
    
Screen('BlendFunction', display.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

Screen('Flip',display.windowPtr);

%% Stimulus Parameters

% p holds stimulus parameters in visual degree.
p.dia = 5;                          % aperture width
p.size = .2;                        % diameter of dots in degrees (.6 Saenz 2003)(.2 Serences 2007)(Ernst .15)
p.nDots = 100;                      % number of dots (?? Serences 2007)(50/field Saenz 2003)(Ernst 200)
p.color = [0 0 0];
p.speed = 5;                        % translational velocity (deg/sec) (4.6 Serences 2007, Ernst)
p.eventspeed = 5;                   % appiditive speed of an event
p.coherence = 1;                    % (Ernst .1)
p.appCentx = 5.9;                   % 5.9 Serences, 11 Saenz
p.appCenty = 2.75;                  % +2.75 Serences, -2.5 Saenz
p.ecc = [-p.appCentx -p.appCentx, p.appCentx p.appCentx;...
    p.appCenty,p.appCenty,p.appCenty,p.appCenty]; % x1,x1,x2,x2;y1,y1,y2,y2 fields 1 and 2
p.lifetime = 12;                    % lifetime in frames (16 Serences 2007) (12 Saenz 2003)
p.direction = [45,135,45,135];      % direction clockwise from up. fields 1 and 2

% dots structure used to hold temporary conversion to pixel locations.
pixecc = angle2pix(display,p.ecc);
dots.ecc(1,:) = [repmat(pixecc(1,1),1,2*p.nDots),repmat(pixecc(1,3),1,2*p.nDots)];
dots.ecc(2,:) = [repmat(pixecc(2,1),1,2*p.nDots),repmat(pixecc(2,3),1,2*p.nDots)];
dots.size = angle2pix(display,p.size);
dots.life = ceil(rand(1,4*p.nDots)*p.lifetime);

% Using centered dots then shifting for each patch so these define
% aperture edges.
left = 0-p.dia/2;
right = 0+p.dia/2;
bottom = 0-p.dia/2;
top = 0+p.dia/2;

% indexing for each patch
L1 = 1:p.nDots;
L2 = (p.nDots+1):2*p.nDots;
R1 = (2*p.nDots+1):3*p.nDots;
R2 = (3*p.nDots+1):4*p.nDots;

%% Timing (secs then frames)
time.pretrial = .2;         % fixation time before cue
time.pretrialframes = round(time.pretrial/display.ifi);
time.cue = 1;              % cue time before dots
time.cueframes = round(time.cue/display.ifi);
time.CSI = .5;              % time between cue and dots
time.CSIframes = round(time.cue/display.ifi);
time.ITI = .3;              % time between end of last trial and start of next
time.ITIframes = round(time.ITI/display.ifi);
time.stim = 3;              % how long for dots
time.stimframes = round(time.stim/display.ifi);
time.resp = .5;             % how long to respond (window to respond in, does not make trial happen faster)
time.respframes = round(time.resp/display.ifi);
time.TT = time.pretrial + time.ITI + time.stim + time.resp;

t = linspace(0,time.stim,time.stimframes+1)';
t = t(1:end-1); % make the right length


%% Event timing
if ~meta.debugging
    cond.event = round(rand(4,cond.nTrials*cond.nBlocks)); % Matrix of zeros and ones for whether an event occured in each surface
else
    cond.event = ones(4,cond.nTrials*cond.nBlocks);  % for debugging!
end

p.tsdsec = .1;              % target temporal wavefore stand. dev. in sec
p.eventbuffer = 3*p.tsdsec; % Event cant occur in the first or last few frames.

% if an event occurs these are the times of the peak:
cond.eventtime = p.eventbuffer + (time.stim - 2*p.eventbuffer)*rand(4,cond.nTrials*cond.nBlocks);

% Speed vector for each frame set to baseline speed
framespeed = p.speed*ones(time.stimframes,1); % Most of the time the dots update with this speed

%% Cue
p.cue.size = .5; % deg width of square
p.cue.color = [0,0,0];
p.cue.ecc = [-1,1];                 % deg ecc from fix
p.cue.loc = [0 0 .25 .15];          % deg

cue.ecc = angle2pix(display,p.cue.ecc);
cue.loc = angle2pix(display,p.cue.loc);

% Cue locations
cueLeftLoc = CenterRectOnPointd(cue.loc, ...
    cue.ecc(1) + display.center(1), 0+display.center(2));
cueRightLoc = CenterRectOnPointd(cue.loc, ...
    cue.ecc(2) + display.center(1), 0+display.center(2));

dimension_array = [25,25,10,5]; % Dimensions of the arrow: head_width, head_height, body_width, body_height

%%
try
    % HideCursor;
    ListenChar(2);
    
    %% Dots
    % Generate the dot starting locations centered on the screen center for
    % all patches
    dots.xy = (rand(2,4*p.nDots)-.5)*p.dia; % (deg) Add eccentricity shift later
    
    %% Start Trial Loop
    for block = 1:cond.nBlocks
        for trial = 1:cond.nTrials
            t = trial + cond.nTrials*(block-1); % count for total trials
            tvec = repmat(framespeed,1,4);
            
            % Set up event
            for field = 1:4 % Need to check if there is an event for every field
                
                
                
                if cond.event(field,trial)
                    % calculate temproal waveform (0-1) with tmean and tsd
                    tmean = cond.eventtime(block,trial);
                    
                    tvec(:,field) = tvec(:,field) + p.eventspeed*exp(-((t-cond.eventtime(field,trial)).^2/(2*p.tsdsec^2)));
                    
                    
                    % put stimulus into larger noise patch (using tx,ty)
                    %                 t1 = cond.eventtime(block,trial);
                    %                 t2 = t1 + g_length - 1;
                    %                 framespeed(:,trial) = framespeed(:,trial) + tvec;                % put in stimulus
                    
                end
            end
            FlushEvents
            
            % Instructions if first trial
            if trial == 1
                DrawFormattedText(display.windowPtr, 'Press Any Key To Begin', 'center', 'center', [0,0,0]);
                Screen('Flip', display.windowPtr);
                KbStrokeWait;
            end
            
            % Pretrial fixation
            drawFix(display);
            vbl = Screen('Flip',display.windowPtr);
            for frame = 1:time.pretrialframes-1 % already used one frame for the fixation
                drawFix(display);
                vbl = Screen('Flip',display.windowPtr,vbl+display.adj);
            end
            
            % Precue
            for frame = 1:time.cueframes
               
                
                
                if cond.blockleftdir(block) == 1 % If Left side is up (1)
                    cuedir(1) = p.direction(1);
                    if cond.blockatt(block) == 1 % If block condition is SAME
                    cuedir(2) = p.direction(1); % same attend direction
                    else
                    cuedir(2) = p.direction(2); % different attend direction
                    end
                else % If Left side is attend down (2)
                    cuedir(1) = p.direction(2);
                    if cond.blockatt(block) == 1 % If block condition is SAME
                    cuedir(2) = p.direction(2);
                    else
                    cuedir(2) = p.direction(1);
                    end
                end
                cuedir = cuedir + 90;
                
                drawFix(display);
                draw_arrow(display.windowPtr, cueLeftLoc, cuedir(1), p.cue.color, dimension_array);
                draw_arrow(display.windowPtr, cueRightLoc, cuedir(2), p.cue.color, dimension_array);

                vbl = Screen('Flip',display.windowPtr,vbl+display.adj);
            end
            % Wait between cue and dots
            for frame = 1:time.CSIframes
                drawFix(display);
                vbl = Screen('Flip',display.windowPtr,vbl+display.adj);
            end
            
            Priority(topPriorityLevel);
            for i = 1:time.stimframes
                
                 % Get the speed for this frame
                dx = zeros(1,4);
                dy = zeros(1,4);
                for f = 1:4
                    curspeed = tvec(i,f)/display.frameRate;
                    
                    % Dot speed - i.e. displacement per frame
                    dx(f) = curspeed*sin(p.direction(f)*pi/180);
                    dy(f) = -curspeed*cos(p.direction(f)*pi/180);
                end
                % Unwrap so they can be added as matrices later
                dx = reshape(repmat(dx,p.nDots,1),[1,4*p.nDots]);
                dy = reshape(repmat(dy,p.nDots,1),[1,4*p.nDots]);
                
                
                x = dots.xy(1,:);
                y = dots.xy(2,:);
                
                % Use the equation of an ellipse to determine which dots fall inside.
                goodDots = x.^2/(p.dia/2)^2 + y.^2/(p.dia/2)^2 < 1;
                
                % Convert to pixels and move to desired eccentricity from screen center
                pixpos.x = angle2pix(display,x) + display.resolution(1)/2  + dots.ecc(1,:);
                pixpos.y = angle2pix(display,y) + display.resolution(2)/2  + dots.ecc(2,:);
                
                % Draw dots within the aperture
                Screen('DrawDots',display.windowPtr,...
                    [pixpos.x(goodDots);pixpos.y(goodDots)], ...
                    dots.size, p.color,[0,0],1);
                
                drawFix(display);
                Screen('Flip',display.windowPtr, vbl+display.adj);
                
                % Update for the next frame
                x = x + dx; % (deg) centered on screen center
                y = y + dy;
                
                %move the dots that are outside the aperture back one aperture
                %width.
                x(x<left) = x(x<left) + p.dia;
                x(x>right) = x(x>right) - p.dia;
                y(y<bottom) = y(y<bottom) + p.dia;
                y(y>top) = y(y>top) - p.dia;
                
                %increment the 'life' of each dot
                dots.life = dots.life+1;
                
                %find the 'dead' dots
                deadDots = mod(dots.life,p.lifetime)==0;
                
                %replace the positions of the dead dots to a random location
                x(deadDots) = (rand(1,sum(deadDots))-.5)*p.dia;
                y(deadDots) = (rand(1,sum(deadDots))-.5)*p.dia;
                
                % move current field dots back to dots structure
                dots.xy(1,:) = x;
                dots.xy(2,:) = y;
                
            end
            Screen('Flip',display.windowPtr);
            Priority(0);
            
            %% Collect response
            
            
            WaitSecs(time.ITI)
            
            %% Save each trial
            % what to save:
            data.block(t,:) = block;
%             data.chosenResp(t,:) = response;      % The rating scale value -2,-1,1,2
%             data.ynresponse(t,:) = ynresponse;     % binary yes no response
%             data.respSide(t,:) = curRespSide;     % Which side the response was given for: 1-Left
%             data.respCorrect(t,:) = correct;      % Binary was the response correct
%             data.correctResponse(t,:) = correctresponse; % Is there a target on the response side
%             data.trialTargetR(t,:) = p.trialTargetR(thistrial);    % Was there a change on the R
%             data.trialTargetL(t,:) = p.trialTargetL(thistrial);
%             data.blockCond(t,:) = p.blockVector(thistrial);           % What is the trial type for this block
%             data.orientLone(t,:) = firstorientation(1);    % 256 RGB for the color bblob
%             data.orientRone(t,:) = firstorientation(2);
%             data.orientLtwo(t,:) = secondorientation(1);
%             data.orientRtwo(t,:) = secondorientation(2);
%             data.RT(t,:) = rt;
%             data.gabortime(t,:) = gabortime;    % [first interval, second interval] (frames)
%             data.gaborloc{t,:} = gaborxy;    % gabor locations. Row - each movie (left then right), column, x,y position. 
%             data.abort(t,:) = trialAbort;   % 0 no abort, 1 = abort
            
        end
    end
    
catch ME
    %Close Screen:
    sca
    ShowCursor;
    ListenChar(1);
    
    rethrow(ME)
end

%% Close Screen:
sca
ShowCursor;
ListenChar(1);
