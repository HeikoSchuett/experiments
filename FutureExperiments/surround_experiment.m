function sampling_struct=surround_experiment(ProbandName,sizeDot,propFilled,lum,time,excentr,Ntrials,crt,use_eyetracker,degSize)
% function sampling_struct=surround_experiment(ProbandName,time,excentr,Ntrials,crt,use_eyetracker)
% this function runs a specific measurement for the parasol experiment
%
% input:
% subject   string for which subject was measured
% time     timing
%          1 = lowpass 1 second
%          2 = maximum flicker (black white every frame)
%          3 = half maximum flicker (1 0 -1 0 1 ....)
%          4 = 15 Hz Sine-Wave Flicker
%          5 = 20 Hz Sine-Wave Flicker
%          6 = 10 Hz Sine-Wave Flicker
% size     size of the dot                        [px]
% excentr  distance from fixation point           [deg]
% Ntrials   total number of trials                        [#]
% crt      optional: run on crt
% use_eyetracker optional: turn on eye_tracker check

addpath('adaptive_sampling')
addpath(genpath('~/Dokumente/ishow'))

rng('shuffle');

if ~exist('sizeDot','var') || isempty(sizeDot)
    sizeDot= 5;
end
if ~exist('Ntrials','var') || isempty(Ntrials)
    Ntrials= 50;
end
if ~exist('use_eyetracker','var') || isempty(use_eyetracker)
    use_eyetracker= 0;
end
if ~exist('crt','var') || isempty(crt)
    crt= 0;
end

if ~exist('degSize','var') || isempty(degSize)
    degSize = [2,2];
end

nWedges = 6;
alpha = 0.6;

%% parse input
assert(length(ProbandName)>=3,'ProbandName has to be at least 3 characters long')
%% setup parameters
% hardware configs
bg_color      = 0.2;
aud_volume    = 0.5;  % audio volume [0,1]
if crt==1
    % CRT
    dist_monitor  = 1650;        % mm
    monitor_px    = [1600,1200]; % px
    monitor_mm    = [400,300];   % mm
elseif crt == 2
    % office
    dist_monitor  = 1000;        % mm
    monitor_px    = [1920,1200]; % px
    monitor_mm    = [520,324];   % mm
else
    % VPIXX
    dist_monitor  = 1350;        % mm
    monitor_px    = [1920,1200]; % px
    monitor_mm    = [484,302];   % mm
end



% experiment times
n_feedback    = 20;   % How many trials per feedback screen

nframes_stim  = 121; % present the whole stimulus for how long?
nframes_iti   = 60; % number of frames between trials-> after response
n_frames_resp = 120;

% saving config
%datadir       = '/home/data/heiko/parasol';
datadir       = '/home/heiko/MATLAB/test';
fileTimestamp = standard_now;
eye_data_path = fullfile(datadir, '/eye_data_files/');

% adaptive sampling
options       = struct;
options.possibleX = exp(linspace(log(.005),log(0.64),50));
options.sigmoidName = 'Weibull';

adaptiveType  = 1;
%% eyeTracker
maxFixDist = 1; %deg




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
if crt == 1
    calib = load('/home/data/calibration/crt_gray/crt_gray2018_03_13_1620.mat');
    clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
    clut(1:10)  = 0;
    clut(clut<0)= 0;
elseif crt == 0
    calib = load('/home/data/calibration/lcd_gray/lcd_gray2018_03_13_1620.mat');
    clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
    clut(1:10)  = 0;
    clut(clut<0)= 0;
end

%% start hardware
if crt == 1
    % CRT-initialization
    win = window('crt_gray', 'bg_color', bg_color, 'clut', clut);
    aud = dpixx_audio_port('volume', aud_volume);
    list_wait = listener_buttonbox;
    list_stop = listener_buttonbox('does_interrupt', true);
elseif crt == 2
    %Debugging at office initialization
    win = window('debug', 'rect', [1920 0 2*1920 1200], 'bg_color', bg_color);
    aud = local_audio_port('volume', aud_volume);
    list_wait = listener_keyboard;
    list_stop = listener_keyboard('does_interrupt', true);
else
    %
    %Datapixx-initialization:
    win = window('lcd_gray', 'bg_color', bg_color, 'clut', clut);
    aud = dpixx_audio_port('volume', aud_volume);
    list_wait = listener_buttonbox;
    list_stop = listener_buttonbox('does_interrupt', true);
end


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
        Eyelink('SetAddress','134.2.202.185');
        
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
        
        eye_used = -1;
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

% calculated values
monitor_deg(1)= atan2(monitor_mm(1),dist_monitor)*180/pi;
monitor_deg(2)= atan2(monitor_mm(2),dist_monitor)*180/pi;
[~,idx] = min(monitor_deg);
pxPerDeg      = monitor_px(idx)./monitor_deg(idx);            % in both direction... should be equal...
pxPerDeg = [pxPerDeg,pxPerDeg];
if numel(degSize)== 1
    sizepx       = 2*floor([degSize,degSize].*pxPerDeg/2)+1;
    degSize         = [degSize,degSize];
else % assume it has 2 elements
    sizepx       = 2*floor(degSize.*pxPerDeg/2)+1;
end

sizeCross  = 5;
size0cross = sizepx+2*sizeCross;

posleft   = [cx,cy]+round([-excentr,0].*pxPerDeg-sizepx/2); %upper left corner
posright  = [cx,cy]+round([excentr,0].*pxPerDeg-sizepx/2); %upper left corner

%% Create textures
% generate fixed grating texture
x= (1:sizeDot)-((sizeDot+1)/2);
r = sqrt(bsxfun(@plus,x.^2,x'.^2));
dot = double(r<=(sizeDot/2));
dot_tex = win.make_texture(dot);

%wedges = image_wedgeTarget(sizepx,6,propFilled,sizeDot,sizeDot+2,sizeDot+5,sizepx(1)/2,0,0,0,1);
antialiasing = 8;
imSize = antialiasing*sizepx;
imageBitmap = zeros(imSize);
x = 1:imSize(2);
y = 1:imSize(1);
x = x-floor(imSize(2)/2);
y = y-floor(imSize(1)/2);

x = x./antialiasing;
y = y./antialiasing;

[xx,yy]=meshgrid(x,y');

r = sqrt(xx.^2+yy.^2);
phi = atan2(yy,xx);
phi2 = phi;
phi2(phi2<0) = phi2(phi2<0)+2*pi;

background = (r>0.1) & (r<=sizepx(1)/2);

wedges = zeros(size(imageBitmap));
wAngles = linspace(-pi,pi,nWedges+1);
anglePerWedge = 4*pi/nWedges*propFilled;
for iw = 1:nWedges
    if (wAngles(iw)+anglePerWedge)<=pi
        wedges = wedges+(background.*(phi>=wAngles(iw)).*(phi<(wAngles(iw)+anglePerWedge))...
            .* (0.5-0.5*cos(2.*pi./anglePerWedge.*(phi-wAngles(iw)))));
    else
        wedges = wedges+(background.*(phi2>=wAngles(iw)).*(phi2<(wAngles(iw)+anglePerWedge))...
            .* (0.5-0.5*cos(2.*pi./anglePerWedge.*(phi2-wAngles(iw)))));
    end
    
end
wedges = imresize(wedges,1/antialiasing,'box');
wedges_tex =  win.make_texture(wedges);

cross_grate = patch_crosshair('size', sizepx, ...
    'length', 40, 'color', .5,'thickness',sizeCross);
cross_tex = win.make_texture(cross_grate.patch);



%% Generate sounds
aud.create_beep('short_low' , 'low' , 5/72, 0.25);
aud.create_beep('short_med' , 'med' , 5/72, 0.25);
aud.create_beep('short_high', 'high', .1, 0.25);
aud.create_beep('long_low'  , 'low' , .4, 0.25);
aud.create_beep('long_high' , 'high', .4, 0.25);


%% run experiment
% initialize sampling_struct for adaptive sampling
sampling_struct = init_adaptive_sampling(0.05,[0.005,0.2],'2AFC',options);
try
    results = [sizeDot,propFilled,lum,time,excentr,NaN,NaN,NaN,NaN]; % c, position, position clicked,RT
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
        c = choose_adaptive_sampling(sampling_struct,adaptiveType);
        if c < sampling_struct.last_stim/4
            c = sampling_struct.last_stim/4;
        elseif c>sampling_struct.last_stim*4
            c = sampling_struct.last_stim*4;
        end
        results(itrial,6) = c;
        redoTrial = true;
        irepeat   = 1;
        krepeat   = 0;
        while redoTrial
            pos    = randi(2);
            results(itrial,7) = pos;
            trialName = [num2str(itrial),num2str(5*krepeat+irepeat)];
            response = trial(c,pos);
            if isnumeric(response)
                redoTrial = true;
            else
                % Save the responses
                [press, rt] = response.get_presses('first');
                redoTrial = false;
                if ~isnan(press)
                    results(itrial,8) = press;
                    results(itrial,9) = rt;
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
            correct=results(:,7)==results(:,8);
            n_correct_total=sum(correct);
            n_total=itrial;
            n_correct_last=sum(correct((itrial-n_feedback+1):itrial));
            % give feedback
            win.pause_trial(list_stop,[num2str(itrial), ' trials passed out of ', num2str(size(results,1)),...
                '\n\n overall ', num2str(100 *n_correct_total/n_total,'%2.1f'), ' % correct'...
                '\n\n ', num2str(n_correct_last), ' of the last ', num2str(n_feedback), ' Trials correct'],[],[],[.5,.5,.5]);
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
    function response = trial(c, pos)
        %function trial(c,iphase, pos)
        %
        % Run one trial
        % frequency dependence removed-> only one frequency, instead phase
        % added
        
        %pause_trial % start the trial only when subject presses a button.
        
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
                win.draw_additive(wedges_tex,alpha*lum,[posleft,posleft+sizepx])
                win.draw_additive(wedges_tex,alpha*lum,[posright,posright+sizepx])
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
                win.draw_additive(wedges_tex,alpha*lum,[posleft,posleft+sizepx])
                win.draw_additive(wedges_tex,alpha*lum,[posright,posright+sizepx])
                if pos == 2 % right is correct
                    win.draw_additive(dot_tex,c.*timecourse(itic),[posright+floor(sizepx/2-sizeDot/2),posright+floor(sizepx/2+sizeDot/2)])
                else
                    win.draw_additive(dot_tex,c.*timecourse(itic),[posleft+floor(sizepx/2-sizeDot/2),posleft+floor(sizepx/2+sizeDot/2)])
                end
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
        if trial_valid
            aud.play('short_high');
            % Response interval
            list_wait.start();
            for itic = 1 : n_frames_resp
                win.draw_additive(wedges_tex,alpha*lum,[posleft,posleft+sizepx])
                win.draw_additive(wedges_tex,alpha*lum,[posright,posright+sizepx])
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
            fname = [ProbandName sprintf('_%d_%.1f_%.1f_',sizeDot,propFilled,lum) fileTimestamp '.mat'];
            save(fullfile(datadir, fname), 'results','eyedata_fname','sampling_struct');
        else
            fname = ['test' sprintf('_%d_%.1f_%.1f_',sizeDot,propFilled,lum) fileTimestamp '.mat'];
            save(fullfile(datadir, fname), 'results','eyedata_fname','sampling_struct');
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
