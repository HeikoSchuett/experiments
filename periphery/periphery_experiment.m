function sampling_struct=periphery_experiment(ProbandName,f0,f1,contrast,size0,size1,excentr0,excentr1,Ntrials,background)
% function periphery_experiment(ProbandName,f0,f1,contrast,size0,size1,excentr0,excentr1,Ntrials,background)
% this function runs a specific measurement for the peripheral sensitivety
% and saves the result as a file starting with ProbandName as a name.
%
% input:
% subject   string for which subject was measured
% f0        frequency of the fixed contrast grating       [cyc/deg]
% f1        frequency of the grating changing contrast    [cyc/deg]
% contrast  contrast of the fixed grating                 [Amplitude/L]
% size0     size of the fixed grating                     [deg]
% size1     size of the variable grating                  [deg]
% excentr0  distance from fixation point (fixed)          [deg]
% excentr1  distance from fixation point (variable)       [deg]
% Ntrials   total number of trials                        [#]
%
% background-> atm. nothing, shall code for the background image to add to
% both images

rng('shuffle');
%% parse input

if ~exist('background','var') || isempty(background)
    iBack = 0;
end

%% setup parameters
% hardware configs
bg_color      = 0.5;
aud_volume    = 0.5;  % audio volume [0,1]
% % CRT
% dist_monitor  = 1650;        % mm
% monitor_px    = [1600,1200]; % px
% monitor_mm    = [400,300];   % mm
% % VPixx
dist_monitor  = 900;        % mm
monitor_px    = [1920,1200]; % px
monitor_mm    = [484,302];   % mm
% office
% dist_monitor  = 1650;        % mm
% monitor_px    = [1920,1200]; % px
% monitor_mm    = [520,324];   % mm

%linearization/calibration file
calib = load('/home/data/calibration/lcd_gray/lcd_gray2014_10_13_1713.mat');



% experiment times
n_feedback    = 20;   % How many trials per feedback screen
nframes_stim  = 50; % present the whole stimulus for how long?
nframes_ramp  = 25; % how many frames for the cosine ramp in and ramp out?
nframes_iti   = 50; % number of frames between trials-> after response
n_frames_resp = 100;

% orientations and phases
theta0        = 0.25;  % horizontal
theta1        = 0.25;  % horizontal
nphase        = 12;
phase         = linspace(0,(1-1/nphase),nphase);


% saving config
%datadir       = '/home/data/heiko/short';
datadir       = '/home/heiko/MATLAB/test';
fileTimestamp = standard_now;
eye_data_path = fullfile(datadir, '/eye_data_files/');

% adaptive sampling
options       = struct;
options.possibleX = exp(linspace(log(.005),log(0.5),50));

adaptiveType  = 1;
%% eyeTracker
maxFixDist = 2; %deg
use_eyetracker = 0;





% calculated values
monitor_deg(1)= atan2(monitor_mm(1),dist_monitor)*180/pi;
monitor_deg(2)= atan2(monitor_mm(2),dist_monitor)*180/pi;
pxPerDeg      = monitor_px./monitor_deg;            % in both direction... should be equal...
if numel(size0)== 1
    size0px       = round([size0,size0].*pxPerDeg);
    size0         = [size0,size0];
else % assume it has 2 elements
    size0px       = round(size0.*pxPerDeg);
end
if numel(size1)== 1
    size1px       = round([size1,size1].*pxPerDeg);
    size1         = [size1,size1];
else % assume it has 2 elements
    size1px       = round(size1.*pxPerDeg);
end

sizeCross  = 5;
size0cross = size0px+2*sizeCross;
size1cross = size0px+2*sizeCross;

pos0left   = round([(monitor_deg(1)/2-excentr0),monitor_deg(2)/2].*pxPerDeg-size0px/2); %upper left corner
pos1left   = round([(monitor_deg(1)/2-excentr1),monitor_deg(2)/2].*pxPerDeg-size1px/2); %upper left corner
pos0right  = round([(monitor_deg(1)/2+excentr0),monitor_deg(2)/2].*pxPerDeg-size0px/2); %upper left corner
pos1right  = round([(monitor_deg(1)/2+excentr1),monitor_deg(2)/2].*pxPerDeg-size1px/2); %upper left corner



% timecourse
% raised cosine in and out
% this has exactly nframes_stim at 1 gain and nframes_ramp at beginning and
% end with gains between 0 and 1
timecourse = [.5*cos(pi*(nframes_ramp-(0:(nframes_ramp-1)))./(nframes_ramp+1))+.5,ones([1,nframes_stim]),.5*cos(pi*(1:nframes_ramp)./(nframes_ramp+1))+.5]';


% calibration
clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
clut(1:10)  = 0;
clut(clut<0)= 0;

%% start hardware

% CRT-initialization
%win = window('crt_gray', 'bg_color', bg_color, 'clut', clut);
%aud = dpixx_audio_port('volume', aud_volume);
%list_wait = listener_buttonbox;
%list_stop = listener_buttonbox('does_interrupt', true);

%
%Datapixx-initialization:
win = window('lcd_gray', 'bg_color', bg_color, 'clut', clut);
aud = dpixx_audio_port('volume', aud_volume);
list_wait = listener_buttonbox;
list_stop = listener_buttonbox('does_interrupt', true);

% %Debugging at office initialization
% win = window('debug', 'rect', [0 0 1920 1200], 'bg_color', bg_color);
% aud = local_audio_port('volume', aud_volume);
% list_wait = listener_keyboard;
% list_stop = listener_keyboard('does_interrupt', true);


%% setup eyetracking data files
session = 1;
eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',subj,'_session_', num2str(session), '.edf'));
while exist(eyedata_fname,'file')
    session = session+1;
    eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',subj,'_session_', num2str(session), '.edf'));
end
% the remote file is stored on MS-DOS, so has filename length restrictions:
edf_file_remote = strcat('s-',subj,'-', num2str(session), '.edf');
%




%% start eyetracker
if use_eyetracker
    % initialise eyetracker
    try
        Eyelink('SetAddress','134.2.202.185')
        
        % Provide Eyelink with details about the graphics environment
        % and perform some initializations. The information is returned
        % in a structure that also contains useful defaults
        % and control codes (e.g. tracker state bit and Eyelink key values).
        el=EyelinkInitDefaults(win.h);
        
        % Initialization of the connection with the Eyelink Gazetracker.
        Eyelink('Initialize','PsychEyelinkDispatchCallback');
        
        Eyelink('Initialize','PsychEyelinkDispatchCallback');
        
        [~, vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
        % make sure that we get gaze data from the Eyelink
        Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
        
        eye_used = -1;
        eyedata_fname = fullfile(datadir, strcat('eyeData',ProbandName,fileTimestamp,'.edf'));
    catch eyelink_error
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        % Close window:
        sca;
        % Restore keyboard output to Matlab:
        ListenChar(0);
        commandwindow;
        disp(eyelink_error);
        disp(eyelink_error.message);
        disp(eyelink_error.stack.line);
    end
end



%% fixation target (According to Giessen)
colorOval = 0;%[0 0 0]; % color of the two circles [R G B]
colorCross = 1;%[255 255 255]; % color of the Cross [R G B]

d1 = 0.6; % diameter of outer circle (degrees)
d2 = 0.2; % diameter of inner circle (degrees)

% position
[cx, cy] = RectCenter(win.rect);

%% Create textures
% generate fixed grating texture
for iphase = 1:length(phase)
    grating0 = grating_sine('size',size0px,'frequency',f0.*size0(1),'phase',phase(iphase),'theta',theta0);
    texGrat0(iphase) = win.make_texture(grating0);
    grating1 = grating_sine('size',size1px,'frequency',f1.*size1(1),'phase',phase(iphase),'theta',theta1);
    texGrat1(iphase) = win.make_texture(grating1);
end

cross_grate0 = patch_crosshair('size', size0px, ...
    'length', 40, 'color', 0,'thickness',sizeCross);
cross_tex0 = win.make_texture(cross_grate0.patch);
cross_grate1 = patch_crosshair('size', size1px, ...
    'length', 40, 'color', 0,'thickness',sizeCross);
cross_tex1 = win.make_texture(cross_grate1.patch);

% win_grate = distrib_tapered_cosine('size', img_size, ...
%                                    'alpha', .3);
window_tex0 = win.make_mask_texture('tapered_cosine', ...
    'size', size0px, ...
    'alpha', .3);
window_tex1 = win.make_mask_texture('tapered_cosine', ...
    'size', size1px, ...
    'alpha', .3);

%% Generate sounds
aud.create_beep('short_low' , 'low' , 5/72, 0.25);
aud.create_beep('short_med' , 'med' , 5/72, 0.25);
aud.create_beep('short_high', 'high', .1, 0.25);
aud.create_beep('long_low'  , 'low' , .4, 0.25);
aud.create_beep('long_high' , 'high', .4, 0.25);


%% run experiment
% initialize sampling_struct for adaptive sampling
sampling_struct = init_adaptive_sampling(contrast,[eps,max(.2,2*contrast)],'equalAsymptote',options);
try
    if ~exist('results','var')
        results = [f0,f1,contrast,size0,size1,excentr0,excentr1,iBack,NaN,NaN,NaN,NaN,NaN];
        results = repmat(results,Ntrials,1);
        firstTrial=1;
    else
        warning('You continue an older experiment. We did not test anything about it')
        firstTrial = find(isnan(results(:,9)),1);
        if isempty(firstTrial)
            error('Experiment already finished')
        end
    end
    pause_trial
    if use_eyetracker
        try
            % open file to record data to
            Eyelink('Openfile', edf_file_remote);
            
            % Calibrate the eye tracker
            EyelinkDoTrackerSetup(el);
            
            %         % do a final check of calibration using driftcorrection
            %         EyelinkDoDriftCorrection(el);
            
        catch eyelink_error
            % Shutdown Eyelink:
            Eyelink('Shutdown');
            % Close window:
            sca;
            % Restore keyboard output to Matlab:
            ListenChar(0);
            commandwindow;
            disp(eyelink_error);
            disp(eyelink_error.message);
            disp(eyelink_error.stack.line);
        end
    end
    for itrial = firstTrial:Ntrials
        c = choose_adaptive_sampling(sampling_struct,adaptiveType);
        results(itrial,12) = c;
        pos    = randi(2);
        iphase = randi(length(phase));
        results(itrial,11) = phase(iphase);
        results(itrial,13) = pos;
        redoTrial = true;
        irepeat   = 1;
        while redoTrial
            trialName = [num2str(itrial),num2str(irepeat)];
            response = trial(c,iphase,pos);
            if isnumeric(response)
                redoTrial = true;
            else
                % Save the responses
                [press, rt] = response.get_presses('first');
                redoTrial = false;
                if ~isnan(press)
                    results(itrial,14) = press;
                    results(itrial,15) = rt;
                else
                    redoTrial = true;
                end
            end
            irepeat   = irepeat+1;
        end
        if mod(itrial,n_feedback)==0
            save_data();
            correct=results(:,13)==results(:,14);
            n_correct_total=sum(correct);
            n_total=itrial;
            n_correct_last=sum(correct((itrial-n_feedback+1):itrial));
            % give feedback
            win.pause_trial(list_stop,[num2str(itrial), ' trials passed out of ', num2str(size(results,1)),...
                '\n\n overall ', num2str(100 *n_correct_total/n_total,'%2.1f'), ' % correct'...
                '\n\n ', num2str(n_correct_last), ' of the last ', num2str(n_feedback), ' Trials correct']);
        end
        % update sampling_struct
        correct = press == pos;
        sampling_struct = update_adaptive_sampling(sampling_struct,[c,correct,1]);
        
    end
catch e
    disp('unexpected stop of experiment!')
    wrap_up();
    rethrow(e)
end

%% function definitions


%% trial
    function response = trial(c,iphase, pos)
        %function trial(c,iphase, pos)
        %
        % Run one trial
        % frequency dependence removed-> only one frequency, instead phase
        % added
        
        pause_trial % start the trial only when subject presses a button.
        
        if use_eyetracker
            % start recording eye position
            Eyelink('StartRecording');
            Eyelink('Message', ['start_trial: ',trialName]);
        end
        trial_valid = true;
        
        % ITI
        for itic = 1 : nframes_iti
            if trial_valid
                if pos == 2 % left is changing
                    win.draw(cross_tex0,1,[pos0left-sizeCross,pos0left+size0cross]);
                    win.draw(cross_tex1,1,[pos1right-sizeCross,pos1right+size1cross]);
                else
                    win.draw(cross_tex0,1,[pos0right-sizeCross,pos0right+size0cross]);
                    win.draw(cross_tex1,1,[pos1left-sizeCross,pos1left+size1cross]);
                end
                draw_Fixation
                win.flip();
                if use_eyetracker
                    if ~fix_check
                        trial_valid = false;
                    end
                end
            end
        end
        
        % stimulus presentations
        %aud.play('short_low');
        for itic = 1 : length(timecourse)
            if trial_valid
                if pos == 2 % right is changing
                    win.draw_additive(texGrat0(iphase), contrast * timecourse(itic),[pos0left,pos0left+size0px]);
                    win.draw_mask(window_tex0,[pos0left,pos0left+size0px]);
                    win.draw_additive(texGrat1(iphase), c * timecourse(itic),[pos1right,pos1right+size1px]);
                    win.draw_mask(window_tex0,[pos1right,pos1right+size1px]);
                    win.draw(cross_tex0,1,[pos0left-sizeCross,pos0left+size0cross]);
                    win.draw(cross_tex1,1,[pos1right-sizeCross,pos1right+size1cross]);
                else
                    win.draw_additive(texGrat0(iphase), contrast * timecourse(itic),[pos0right,pos0right+size0px]);
                    win.draw_mask(window_tex0,[pos0right,pos0right+size0px]);
                    win.draw_additive(texGrat1(iphase), c * timecourse(itic),[pos1left,pos1left+size1px]);
                    win.draw_mask(window_tex0,[pos1left,pos1left+size1px]);
                    win.draw(cross_tex0,1,[pos0right-sizeCross,pos0right+size0cross]);
                    win.draw(cross_tex1,1,[pos1left-sizeCross,pos1left+size1cross]);
                end
                draw_Fixation
                win.flip();
                if use_eyetracker
                    if ~fix_check
                        trial_valid = false;
                    end
                end
            end
        end
        if trial_valid
            aud.play('short_high');
            % Response interval
            list_wait.start();
            for itic = 1 : n_frames_resp
                if pos == 1
                    win.draw(cross_tex0,1,[pos0left-sizeCross,pos0left+size0cross]);
                    win.draw(cross_tex1,1,[pos1right-sizeCross,pos1right+size1cross]);
                else
                    win.draw(cross_tex0,1,[pos0right-sizeCross,pos0right+size0cross]);
                    win.draw(cross_tex1,1,[pos1left-sizeCross,pos1left+size1cross]);
                end
                win.flip();
            end
            response = list_wait.stop();
        else
            response = NaN;
            Eyelink('Message', 'trial invalid');
        end
        
        if use_eyetracker
            Eyelink('Message', 'end_trial');
            Eyelink('StopRecording');
        end
    end

%% pause (show only fixation cross)
    function pause_trial
        list_stop.start();
        draw_Fixation;
        win.flip;
        while true
            try
                list_stop.check();
            catch e
                if ~strcmp(e.identifier, 'iShow:ResponseInterrupt')
                    rethrow(e)
                end
                break
            end
        end
    end
%% Fixation cross drawing (from Paper coppied)
    function draw_Fixation
        Screen('FillOval', win.h, colorOval, [cx-d1/2 * pxPerDeg(1), cy-d1/2 * pxPerDeg(1), cx+d1/2 * pxPerDeg(1), cy+d1/2 * pxPerDeg(1)], d1 * pxPerDeg(1));
        Screen('DrawLine', win.h, colorCross, cx-d1/2 * pxPerDeg(1), cy, cx+d1/2 * pxPerDeg(1), cy, d2 * pxPerDeg(1));
        Screen('DrawLine', win.h, colorCross, cx, cy-d1/2 * pxPerDeg(1), cx, cy+d1/2 * pxPerDeg(1), d2 * pxPerDeg(1));
        Screen('FillOval', win.h, colorOval, [cx-d2/2 * pxPerDeg(1), cy-d2/2 * pxPerDeg(1), cx+d2/2 * pxPerDeg(1), cy+d2/2 * pxPerDeg(1)], d2 * pxPerDeg(1));
    end

%% Fixation check
    function result = fix_check
        result = true;
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            % if we do, get current gaze position from sample
            x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
            y = evt.gy(eye_used+1);
            
            fixCheckX = abs(x-cx)./pxPerDeg(1);
            fixCheckY = abs(y-cy)./pxPerDeg(2);
            
            % do we have valid data and is the pupil visible?
            if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                if sqrt(fixCheckX.^2 + fixCheckY.^2) > maxFixDist
                    result = false;
                end
            end
        end
    end
%% save
    function save_data
        if ~exist(datadir,'dir')
            mkdir(datadir);
        end
        if exist('ProbandName','var')
            fname = [ProbandName fileTimestamp];
            save(fullfile(datadir, fname), 'results');
        else
            fname = ['test' fileTimestamp];
            save(fullfile(datadir, fname), 'results');
        end
        fprintf('Data saved in: %s\n', fname);
    end
%% wrap up

    function wrap_up
        save_data
        try
            if use_eyetracker
                Eyelink('Command', 'clear_screen 0')
                Eyelink('CloseFile');
                % download data file
                try
                    fprintf('Receiving data file ''%s''\n', edf_file_remote );
                    status=Eyelink('ReceiveFile',[],eye_data_path,1);
                    pause(1);
                    
                    movefile([eye_data_path, edf_file_remote], eyedata_fname);
                    if status > 0
                        fprintf('ReceiveFile status %d\n', status);
                    end
                    if 2==exist(edf_file_remote, 'file')
                        fprintf('Data file ''%s'' can be found in ''%s''\n',...
                        edf_file_remote, eye_data_path );
                    end
                catch
                    fprintf('Problem receiving data file ''%s''\n', edf_file_remote );
                end
                % Shutdown Eyelink:
                Eyelink('Shutdown');
            end
            win.clear;
            win.delete;
            clear win
        catch
            Screen('CloseAll')
        end
        %         plot_experiment(results);
        %         disp('   contrast, Ncorrect,    N');
        %         disp(accumulate_trials(results(:,9:11)));
    end
end