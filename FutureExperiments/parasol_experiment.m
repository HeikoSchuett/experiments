function sampling_struct=parasol_experiment(ProbandName,f,time,type,degSize,excentr,Ntrials,crt,use_eyetracker)
% function parasol_experiment(ProbandName,f,time,size0,excentr0,Ntrials,crt,use_eyetracker)
% this function runs a specific measurement for the parasol experiment
%
% input:
% subject   string for which subject was measured
% f        frequency of the grating       [cyc/deg]
% time     timing 
%          1 = lowpass 1 second
%          2 = maximum flicker (black white every frame)
%          3 = half maximum flicker (1 0 -1 0 1 ....)
%          4 = 15 Hz Sine-Wave Flicker
%          5 = 20 Hz Sine-Wave Flicker
%          6 = 10 Hz Sine-Wave Flicker
% type     1 = detection
%          2 = orientation discrimination
%          3 = orientation discrimination near horizontal
% size     size of the  grating                     [deg]
% excentr  distance from fixation point           [deg]
% Ntrials   total number of trials                        [#]
% crt      optional: run on crt
% use_eyetracker optional: turn on eye_tracker check
% 
% to convert to database format use:
% convert_database('data','output.csv',{'frequency','temporal condition','detection condition','sizeX','sizeY','excentricity','test contrast','phase','true position','response','RT'})
    
addpath('adaptive_sampling')
addpath(genpath('~/Dokumente/ishow'))

rng('shuffle');

if ~exist('use_eyetracker','var') || isempty(use_eyetracker)
    use_eyetracker= 0;
end
if ~exist('crt','var') || isempty(crt)
    crt= 0;
end

%% parse input
assert(length(ProbandName)>=3,'ProbandName has to be at least 3 characters long')
%% setup parameters
% hardware configs
bg_color      = 0.5;
aud_volume    = 0.5;  % audio volume [0,1]
if crt
    % % CRT
    dist_monitor  = 1650;        % mm
    monitor_px    = [1600,1200]; % px
    monitor_mm    = [400,300];   % mm
else
    % VPIXX
    dist_monitor  = 1350;        % mm
    monitor_px    = [1920,1080]; % px
    monitor_mm    = [534,301];   % mm
end
% office
% dist_monitor  = 1000;        % mm
% monitor_px    = [1920,1200]; % px
% monitor_mm    = [520,324];   % mm




% experiment times
n_feedback    = 20;   % How many trials per feedback screen

nframes_stim  = 121; % present the whole stimulus for how long?
nframes_iti   = 60; % number of frames between trials-> after response
n_frames_resp = 120;

% orientations and phases
theta0        = pi/2;  % horizontal
nphase        = 12;
phase         = linspace(0,(1-1/nphase),nphase);


% saving config
%datadir       = '/home/data/heiko/parasol';
%datadir       = '/home/heiko/MATLAB/test';
datadir       = '/home/wiebke/results';
fileTimestamp = standard_now;
eye_data_path = fullfile(datadir, '/eye_data_files/');

% adaptive sampling
options       = struct;
options.possibleX = exp(linspace(log(.005),log(0.64),50));
options.sigmoidName = 'Weibull';

adaptiveType  = 1;
%% eyeTracker
maxFixDist = 1; %deg
eye_used = -1;

% calculated values
monitor_deg(1)= atan2(monitor_mm(1),dist_monitor)*180/pi;
monitor_deg(2)= atan2(monitor_mm(2),dist_monitor)*180/pi;
[~,idx] = min(monitor_deg);
pxPerDeg      = monitor_px(idx)./monitor_deg(idx);            % in both direction... should be equal...
pxPerDeg = [pxPerDeg,pxPerDeg];
if numel(degSize)== 1
    sizepx       = round([degSize,degSize].*pxPerDeg);
    degSize         = [degSize,degSize];
else % assume it has 2 elements
    sizepx       = round(degSize.*pxPerDeg);
end

sizeCross  = 5;
size0cross = sizepx+2*sizeCross;

%% fixation target (According to Giessen)
colorOval = 0;%[0 0 0]; % color of the two circles [R G B]
colorCross = 1;%[255 255 255]; % color of the Cross [R G B]
% timecourse
% raised cosine in and out
timecourse = linspace(-pi,pi,nframes_stim)';
timecourse = 0.5+0.5*cos(timecourse);
switch time
    case 1 % timecourse is done
    case 2
        timecourse = timecourse.*[repmat([1;-1],floor(nframes_stim/2),1);1];
    case 3
        timecourse = timecourse.*[repmat([1;0;-1;0],floor(nframes_stim/4),1);1];
    case 4
        timecourse = timecourse.*cos(15*linspace(-pi,pi,nframes_stim)');
    case 5
        timecourse = timecourse.*cos(20*linspace(-pi,pi,nframes_stim)');
    case 6
        timecourse = timecourse.*cos(10*linspace(-pi,pi,nframes_stim)');
end

% calibration
%linearization/calibration file
if crt
    calib = load('/home/data/calibration/crt_gray/crt_gray2018_03_13_1620.mat');
else
    calib = load('/home/data/calibration/lcd_gray/lcd_gray2019_04_29_1105.mat');
end
clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
clut(1:10)  = 0;
clut(clut<0)= 0;

%% start hardware
if crt
    % CRT-initialization
    win = window('crt_gray', 'bg_color', bg_color, 'clut', clut);
    aud = dpixx_audio_port('volume', aud_volume);
    list_wait = listener_buttonbox;
    list_stop = listener_buttonbox('does_interrupt', true);
else
    %
    %Datapixx-initialization:
    win = window('lcd_gray', 'bg_color', bg_color, 'clut', clut);
    aud = local_audio_port('volume', aud_volume);
    list_wait = listener_buttonbox;
    list_stop = listener_buttonbox('does_interrupt', true);
end
% %Debugging at office initialization
% win = window('debug', 'rect', [1920 0 2*1920 1200], 'bg_color', bg_color);
% aud = local_audio_port('volume', aud_volume);
% list_wait = listener_keyboard;
% list_stop = listener_keyboard('does_interrupt', true);


%% setup eyetracking data files
session = 1;
eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',ProbandName,'_session_', num2str(session), '.edf'));
while exist(eyedata_fname,'file')
    session = session+1;
    eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',ProbandName,'_session_', num2str(session), '.edf'));;
end
% the remote file is stored on MS-DOS, so has filename length restrictions:
edf_file_remote = strcat('s_',ProbandName(1:3),'_', num2str(session), '.edf');
%




%% start eyetracker
if use_eyetracker
    % initialise eyetracker
    try
        Eyelink('SetAddress','192.168.178.101');
        
        % Provide Eyelink with details about the graphics environment
        % and perform some initializations. The information is returned
        % in a structure that also contains useful defaults
        % and control codes (e.g. tracker state bit and Eyelink key values).
        el=EyelinkInitDefaults(win.h);
        
        % Initialization of the connection with the Eyelink Gazetracker.
        Eyelink('Initialize','PsychEyelinkDispatchCallback');
        
        
        [~, vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
        % make sure that we get gaze data from the Eyelink
        Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
       
    catch eyelink_error
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        % Close window:
        sca;
        % Restore keyboard output to Matlab:
        ListenChar(0);
        commandwindow;
        disp(eyelink_error);
        %disp(eyelink_error.message);
        %disp(eyelink_error.stack.line);
    end
end



%% fixation target (According to Giessen)
colorOval = 0;%[0 0 0]; % color of the two circles [R G B]
colorCross = 1;%[255 255 255]; % color of the Cross [R G B]

d1 = 0.6; % diameter of outer circle (degrees)
d2 = 0.2; % diameter of inner circle (degrees)

% position
[cx, cy] = RectCenter(win.rect);

sizeCross  = 5;
size0cross = sizepx+2*sizeCross;

posleft   = [cx,cy]+round([-excentr,0].*pxPerDeg-sizepx/2); %upper left corner
posright  = [cx,cy]+round([excentr,0].*pxPerDeg-sizepx/2); %upper left corner
posleft   = round([(monitor_deg(1)/2-excentr),monitor_deg(2)/2].*pxPerDeg-sizepx/2); %upper left corner
posright  = round([(monitor_deg(1)/2+excentr),monitor_deg(2)/2].*pxPerDeg-sizepx/2); %upper left corner




%% Create textures
% generate fixed grating texture
for iphase = 1:length(phase)
    if (type == 1) || (type == 2)
        grating = grating_sine('size',sizepx,'frequency',f.*degSize(1),'phase',phase(iphase),'theta',theta0);
    else
        grating = grating_sine('size',sizepx,'frequency',f.*degSize(1),'phase',phase(iphase),'theta',theta0+pi/8);
    end
    if type == 1
        gratingWrong = zeros(sizepx);
    elseif type == 2
        gratingWrong = grating_sine('size',sizepx,'frequency',f.*degSize(1),'phase',phase(iphase),'theta',theta0+pi/2);
    else
        gratingWrong = grating_sine('size',sizepx,'frequency',f.*degSize(1),'phase',phase(iphase),'theta',theta0-pi/8);
    end
    texGrat(iphase,1) = win.make_texture(grating);
    texGrat(iphase,2) = win.make_texture(-grating);
    texGratWrong(iphase,1) = win.make_texture(gratingWrong);
    texGratWrong(iphase,2) = win.make_texture(-gratingWrong);
end

cross_grate = patch_crosshair('size', sizepx, ...
    'length', 40, 'color', 0,'thickness',sizeCross);
cross_tex = win.make_texture(cross_grate.patch);

% win_grate = distrib_tapered_cosine('size', img_size, ...
%                                    'alpha', .3);
window_tex = win.make_mask_texture('tapered_cosine', ...
    'size', sizepx, ...
    'alpha', 0.2);

%% Generate sounds
aud.create_beep('short_low' , 'low' , 5/72, 0.25);
aud.create_beep('short_med' , 'med' , 5/72, 0.25);
aud.create_beep('short_high', 'high', .1, 0.25);
aud.create_beep('long_low'  , 'low' , .4, 0.25);
aud.create_beep('long_high' , 'high', .4, 0.25);


%% run experiment
% initialize sampling_struct for adaptive sampling
sampling_struct = init_adaptive_sampling(0.1,[0.0005,0.75],'2AFC',options);
try
    results = [f,time,type,degSize,excentr,NaN,NaN,NaN,NaN,NaN,NaN]; % c, phase, position, position clicked,RT
    results = repmat(results,Ntrials,1);
    firstTrial=1;
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
        [~, ~, keyCode] = KbCheck();
        if keyCode(KbName('q'))
            error('Experiment ended by q press')
        end
        c = choose_adaptive_sampling(sampling_struct,adaptiveType);
        if c < sampling_struct.last_stim/4
            c = sampling_struct.last_stim/4;
        elseif c>sampling_struct.last_stim*4
            c = sampling_struct.last_stim*4;
        end
        results(itrial,12) = now;
        results(itrial,7) = c;
        iphase = randi(length(phase));
        results(itrial,8) = phase(iphase);
        redoTrial = true;
        irepeat   = 1;
        krepeat   = 0;
        while redoTrial
            pos    = randi(2);
            results(itrial,9) = pos;
            trialName = [num2str(itrial),num2str(5*krepeat+irepeat)];
            response = trial(c.*bg_color,iphase,pos);
            if isnumeric(response)
                redoTrial = true;
            else
                % Save the responses
                [press, rt] = response.get_presses('first');
                redoTrial = false;
                if ~isnan(press)
                    results(itrial,10) = press;
                    results(itrial,11) = rt;
                else
                    redoTrial = true;
                end
            end
            if irepeat>5
                % Calibrate the eye tracker
                if use_eyetracker
                    EyelinkDoTrackerSetup(el);
                    EyelinkDoDriftCorrection(el);
                end
                irepeat = 0;
                krepeat = krepeat+1;
            end
            irepeat   = irepeat+1;
        end
        if mod(itrial,n_feedback)==0
            save_data();
            correct=results(:,9)==results(:,10);
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

wrap_up();


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
        vbl=GetSecs;
        for itic = 1 : nframes_iti
            if trial_valid
                win.draw(cross_tex,1,[posleft-sizeCross,posleft+size0cross]);
                win.draw(cross_tex,1,[posright-sizeCross,posright+size0cross]);
                draw_Fixation
                vbl = Screen('Flip', win.h,vbl+0.5/win.framerate);
                if use_eyetracker
                    if ~fix_check
                        trial_valid = false;
                    end
                end
            end
        end
        
        % stimulus presentations
        aud.play('short_low');
        for itic = 1 : length(timecourse)
            if trial_valid
                if pos == 2 % right is correct
                    if c * timecourse(itic)>0
                        win.draw_additive(texGratWrong(iphase,1), c * timecourse(itic),[posleft,posleft+sizepx]);
                        win.draw_additive(texGrat(iphase,1), c * timecourse(itic),[posright,posright+sizepx]);
                    else
                        win.draw_additive(texGratWrong(iphase,2), -c * timecourse(itic),[posleft,posleft+sizepx]);
                        win.draw_additive(texGrat(iphase,2), -c * timecourse(itic),[posright,posright+sizepx]);
                    end
                    win.draw_mask(window_tex,[posleft,posleft+sizepx]);
                    win.draw_mask(window_tex,[posright,posright+sizepx]);
                    win.draw(cross_tex,1,[posleft-sizeCross,posleft+size0cross]);
                    win.draw(cross_tex,1,[posright-sizeCross,posright+size0cross]);
                else
                    if c * timecourse(itic)>0
                        win.draw_additive(texGratWrong(iphase,1), c * timecourse(itic),[posright,posright+sizepx]);
                        win.draw_additive(texGrat(iphase,1), c * timecourse(itic),[posleft,posleft+sizepx]);
                    else
                        win.draw_additive(texGratWrong(iphase,2), -c * timecourse(itic),[posright,posright+sizepx]);
                        win.draw_additive(texGrat(iphase,2), -c * timecourse(itic),[posleft,posleft+sizepx]);
                    end
                    win.draw_mask(window_tex,[posright,posright+sizepx]);
                    win.draw_mask(window_tex,[posleft,posleft+sizepx]);
                    win.draw(cross_tex,1,[posright-sizeCross,posright+size0cross]);
                    win.draw(cross_tex,1,[posleft-sizeCross,posleft+size0cross]);
                end
                draw_Fixation
                vbl = Screen('Flip', win.h,vbl+0.5/win.framerate);
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
                win.draw(cross_tex,1,[posleft-sizeCross,posleft+size0cross]);
                win.draw(cross_tex,1,[posright-sizeCross,posright+size0cross]);
                vbl = Screen('Flip', win.h,vbl+0.5/win.framerate);
            end
            response = list_wait.stop();
        else
            response = NaN;
            if use_eyetracker
                Eyelink('Message', 'trial invalid');
            end
        end
        
        if use_eyetracker
            Eyelink('Message', 'end_trial');
            Eyelink('StopRecording');
            if ~trial_valid
                WaitSecs(0.5);
            end
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
        %Screen('DrawLine', win.h, colorCross, cx-d1/2 * pxPerDeg(1), cy, cx+d1/2 * pxPerDeg(1), cy, d2 * pxPerDeg(1));
        Screen('FillRect', win.h, colorCross, [cx-d1/2 * pxPerDeg(1), cy-(d2/2 * pxPerDeg(1)), cx+d1/2 * pxPerDeg(1), cy+(d2/2 * pxPerDeg(1))]);
        %Screen('DrawLine', win.h, colorCross, cx, cy-d1/2 * pxPerDeg(1), cx, cy+d1/2 * pxPerDeg(1), d2 * pxPerDeg(1));
        Screen('FillRect', win.h, colorCross, [cx-(d2/2 * pxPerDeg(1)), cy-d1/2 * pxPerDeg(1), cx+(d2/2 * pxPerDeg(1)), cy+d1/2 * pxPerDeg(1)]);
        Screen('FillOval', win.h, colorOval, [cx-d2/2 * pxPerDeg(1), cy-d2/2 * pxPerDeg(1), cx+d2/2 * pxPerDeg(1), cy+d2/2 * pxPerDeg(1)], d2 * pxPerDeg(1));
    end

%% Fixation check
    function result = fix_check
        result = true;
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            
            if eye_used == -1
                % if we do, get current gaze position from sample
                x = mean(evt.gx); % +1 as we're accessing MATLAB array
                y = mean(evt.gy);
                
                fixCheckX = abs(x-cx)./pxPerDeg(1);
                fixCheckY = abs(y-cy)./pxPerDeg(2);
                
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && all(evt.pa>0)
                    if sqrt(fixCheckX.^2 + fixCheckY.^2) > maxFixDist
                        result = false;
                    end
                end
            else
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
    end
%% save
    function save_data
        if ~exist(datadir,'dir')
            mkdir(datadir);
        end
        if exist('ProbandName','var')
            fname = [ProbandName sprintf('_%d_%d_%.1f_',type,time,f) fileTimestamp '.mat'];
            save(fullfile(datadir, fname), 'results','eyedata_fname','sampling_struct');
        else
            fname = ['test' sprintf('_%d_%d_%.1f_',type,time,f) fileTimestamp '.mat'];
            save(fullfile(datadir, fname), 'results','eyedata_fname','sampling_struct');
        end
        fprintf('Data saved in: %s\n', fname);
    end
%% wrap up

    function wrap_up
        save_data
        try
            PsychPortAudio('Stop', aud.h);
            PsychPortAudio('Close', aud.h);
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
            ListenChar(0);
        end
        %         plot_experiment(results);
        %         disp('   contrast, Ncorrect,    N');
        %         disp(accumulate_trials(results(:,9:11)));
    end
end
