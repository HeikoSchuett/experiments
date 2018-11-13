function sampling_struct=witch_experiment(ProbandName,f,typeTemporal,distance,Ntrials,n_block)
% function sampling_struct=witch_experiment(ProbandName,f,typeTemporal,distance,Ntrials,n_block)
% this functions tries to closely reproduce the conditions in Baldwin &
% Meese for the foveal presentation to allow calibration of the Retina
% model
%
% input:
% ProbandName  string for which subject was measured
% f            frequency of the grating(s)                 [cyc/deg]
% typeTemporal temporal condition (1= rectangle, 2= lowpass, 3= highpass,4=100ms,5=1.5 second raised Cosine)
% distance     distance of the monitor                     [mm]
% Ntrials      total number of trials                      [#]
% n_block      number of trials with fixed stimulus intensity (default=1)
%
% please check that framerate = 120hz!

bw = [1.6,25/180*pi];

addpath('adaptive_sampling')
rng('shuffle');
%% parse input

%% setup parameters
% hardware configs
bg_color      = 0.5;
aud_volume    = 0.5;  % audio volume [0,1]
% % CRT
% dist_monitor  = 1650;        % mm
monitor_px    = [1280,960]; % px
monitor_mm    = [400,300];   % mm
framerate = 120; %hz
% % VPixx
%dist_monitor  = 900;        % mm
%monitor_px    = [1920,1200]; % px
%monitor_mm    = [484,302];   % mm
% office
% dist_monitor  = 1650;        % mm
% monitor_px    = [1920,1200]; % px
% monitor_mm    = [520,324];   % mm

dist_monitor = distance; %mm!

n_feedback    = 20;   % How many trials per feedback screen

if ~exist('n_block','var') || isempty(n_block)
    n_block = 1;
end




%linearization/calibration file
%calib = load('/home/data/calibration/crt_gray/crt_gray2016_04_20_1402_1280x960_120.mat');
calib = load('/home/data/calibration/crt_gray/crt_gray2018_09_19_1909.mat');



% experiment times
%nframes_stim  = 50; % present the whole stimulus for how long?
%nframes_ramp  = 25; % how many frames for the cosine ramp in and ramp out?
nframes_iti   = 60; % number of frames between trials-> after response
nframes_isi   = 24;
n_frames_resp = 100;

% timecourse
switch typeTemporal
    case 1 %300ms
        timecourse = [zeros(ceil(.15*framerate),1);ones(ceil(.3*framerate),1);zeros(ceil(.15*framerate),1)];
    case 2 %300ms lowpassfiltered
        timecourse = [zeros(ceil(.15*framerate),1);ones(ceil(.3*framerate),1);zeros(ceil(.15*framerate),1)];
        %timefilter = window(@hann,ceil(.25*framerate)+1);
        NFilter    = ceil(.3*framerate)+1;
        timefilter = 0.5-0.5*cos(2*pi*(0:(NFilter-1))/(NFilter-1));
        timefilter = timefilter./sum(timefilter);
        timecourse = conv(timecourse,timefilter,'same');
    case 3 %300ms highpass filtered
        timecourse = [zeros(ceil(.15*framerate),1);ones(ceil(.3*framerate),1);zeros(ceil(.15*framerate),1)];
        %timefilter = window(@hann,ceil(.25*framerate)+1);
        NFilter    = ceil(.3*framerate)+1;
        timefilter = 0.5-0.5*cos(2*pi*(0:(NFilter-1))/(NFilter-1));
        timefilter = timefilter./sum(timefilter);
        timecourseL = conv(timecourse,timefilter,'same');
        timecourse = timecourse-timecourseL;
    case 4 %witch hat csf
        timecourse = ones(ceil(.1*framerate),1);
    case 5 % 1.5 second raised cosine to copy Felix
        timecourse = (1:1.5*framerate)-1;
        timecourse = 0.5+0.5*sin(2*pi*timecourse./(1.5*framerate-1)-pi/2);
    case 6 %300ms lowpassfiltered + 3 cycles
        timecourse = [zeros(ceil(.15*framerate),1);ones(ceil(.3*framerate),1);zeros(ceil(.15*framerate),1)];
        %timefilter = window(@hann,ceil(.25*framerate)+1);
        NFilter    = ceil(.3*framerate)+1;
        timefilter = 0.5-0.5*cos(2*pi*(0:(NFilter-1))/(NFilter-1));
        timefilter = timefilter./sum(timefilter);
        timecourse = conv(timecourse,timefilter,'same');
        timecourse = timecourse.*-cos(3*2*pi*(0:(length(timecourse)-1))/(length(timecourse)-1))';
    case 7 %300ms lowpassfiltered + 6 cycles
        timecourse = [zeros(ceil(.15*framerate),1);ones(ceil(.3*framerate),1);zeros(ceil(.15*framerate),1)];
        %timefilter = window(@hann,ceil(.25*framerate)+1);
        NFilter    = ceil(.3*framerate)+1;
        timefilter = 0.5-0.5*cos(2*pi*(0:(NFilter-1))/(NFilter-1));
        timefilter = timefilter./sum(timefilter);
        timecourse = conv(timecourse,timefilter,'same');
        timecourse = timecourse.*cos(6*2*pi*(0:(length(timecourse)-1))/(length(timecourse)-1))';
end



% saving config
datadir = '/home/data/heiko/witch';
if ~exist(datadir,'dir')
    mkdir(datadir)
end
%datadir       = '/home/heiko/MATLAB/test';
fileTimestamp = standard_now;

% adaptive sampling
options       = struct;
options.possibleX = exp(linspace(log(.0005),log(.512),51));
options.sigmoidName = 'Weibull';


adaptiveType  = 2;



% calculated values
monitor_deg(1)= 2*atan2(monitor_mm(1)/2,dist_monitor)*180/pi;
monitor_deg(2)= 2*atan2(monitor_mm(2)/2,dist_monitor)*180/pi;
pxPerDeg      = monitor_px./monitor_deg;            % in both direction... should be equal...

size0 = [5,5];

if numel(size0)== 1
    size0px       = round([size0,size0].*pxPerDeg);
    size0         = [size0,size0];
else % assume it has 2 elements
    size0px       = round(size0.*pxPerDeg);
end

size1 = size0;
if numel(size1)== 1
    size1px       = round([size1,size1].*pxPerDeg);
    size1         = [size1,size1];
else % assume it has 2 elements
    size1px       = round(size1.*pxPerDeg);
end

sizeCross  = 3;
size0cross = size0px+2*sizeCross;


pos0  = round([(monitor_deg(1)/2),monitor_deg(2)/2].*pxPerDeg-size0px/2); %upper left corner


% calibration
clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
clut(1:10)  = 0;
clut(clut<0)= 0;

%% start hardware

% CRT-initialization
win = window('crt_gray', 'bg_color', bg_color, 'clut', clut);
aud = dpixx_audio_port('volume', aud_volume);
list_wait = listener_buttonbox;
list_stop = listener_buttonbox('does_interrupt', true);

%
% %Datapixx-initialization:
% win = window('lcd_gray', 'bg_color', bg_color, 'clut', clut);
% aud = dpixx_audio_port('volume', aud_volume);
% list_wait = listener_buttonbox;
% list_stop = listener_buttonbox('does_interrupt', true);

% %Debugging at office initialization
% win = window('debug', 'rect', [0 0 1000 1000], 'bg_color', bg_color);
% aud = local_audio_port('volume', aud_volume);
% list_wait = listener_keyboard;
% list_stop = listener_keyboard('does_interrupt', true);




%% fixation target (According to Giessen)
colorOval = 0;%[0 0 0]; % color of the two circles [R G B]
colorCross = 1;%[255 255 255]; % color of the Cross [R G B]

d1 = 0.6; % diameter of outer circle (degrees)
d2 = 0.2; % diameter of inner circle (degrees)

% position
[cx, cy] = RectCenter(win.rect);

%% Create textures
% generate fixed grating texture
grating0 = image_logGaborOrth(size0px,size0,f,1,size0/2,pi/2,1,bw,0)-1;
texGrat0 = win.make_texture(grating0);

cross_grate0 = patch_crosshair('size', size0px, ...
    'length', 20, 'color', 0,'thickness',sizeCross);
cross_tex0 = win.make_texture(cross_grate0.patch);


%% Generate sounds
aud.create_beep('short_low' , 'low' , length(timecourse)/framerate, 0.25);
aud.create_beep('short_med' , 'med' , length(timecourse)/framerate, 0.25);
aud.create_beep('short_high', 'high', length(timecourse)/framerate, 0.25);
aud.create_beep('long_low'  , 'low' , nframes_iti/2/framerate, 0.25);
aud.create_beep('long_high' , 'high', nframes_iti/2/framerate, 0.25);


%% run experiment
% initialize sampling_struct for adaptive sampling
sampling_struct = init_adaptive_sampling(0.1,options.possibleX,'2AFC',options);
sampling_struct.distance = distance;
sampling_struct.typeTemporal = typeTemporal;
try
    if ~exist('results','var')
        results = [f,size0,NaN,NaN,NaN,NaN];
        results = repmat(results,Ntrials,1);
        firstTrial=1;
    else
        warning('You continue an older experiment. We did not test anything about it')
        firstTrial = find(isnan(results(:,7)),1);
        if isempty(firstTrial)
            error('Experiment already finished')
        end
    end
    blockdata = [];
    pause_trial
    for itrial = firstTrial:Ntrials
        if n_block == 1 ||  mod(itrial,n_block)==1 % update contrast only every Nblockth trial
            c = choose_adaptive_sampling(sampling_struct,adaptiveType);
        end
        results(itrial,4) = c;
        pos    = randi(2);
        results(itrial,5) = pos;
        response = trial(c,pos);
        % Save the responses
        [press, rt] = response.get_presses('first');
        if ~isnan(press)
            results(itrial,6) = press;
            results(itrial,7) = rt;
        else
            results(itrial,6) = NaN;
            results(itrial,7) = NaN;
        end
        % update sampling_struct
        correct = press == pos;
        blockdata = [blockdata;[c,correct,1]];
        if mod(itrial,n_block)==0
            sampling_struct = update_adaptive_sampling(sampling_struct,blockdata);
            blockdata = [];
        end
        if mod(itrial,n_feedback)==0
            save_data();
            correct=results(:,5)==results(:,6);
            n_correct_total=sum(correct);
            n_total=itrial;
            n_correct_last=sum(correct((itrial-n_feedback+1):itrial));
            % give feedback
            win.pause_trial(list_stop,[num2str(itrial), ' trials passed out of ', num2str(size(results,1)),...
                '\n\n overall ', num2str(100 *n_correct_total/n_total,'%2.1f'), ' % correct'...
                '\n\n ', num2str(n_correct_last), ' of the last ', num2str(n_feedback), ' Trials correct']);
        end
        
    end
    save_data();
    correct=results(:,5)==results(:,6);
    n_correct_total=sum(correct);
    n_total=itrial;
    win.pause_trial(list_stop,[num2str(itrial), ' trials passed out of ', num2str(size(results,1)),...
        '\n\n overall ', num2str(100 *n_correct_total/n_total,'%2.1f'), ' % correct']);
catch e
    disp('unexpected stop of experiment!')
    rethrow(e)
    wrap_up();
end

%% function definitions


%% trial
    function response = trial(c, pos)
        %function trial(c,iphase, pos)
        %
        % Run one trial
        % frequency dependence removed-> only one frequency, instead phase
        % added
        
        %pause_trial % start the trial only when subject presses a button.
        
        trial_valid = true;
        
        % ITI
        for itic = 1 : nframes_iti
            if trial_valid
                win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                win.flip();
            end
        end
        
        % stimulus presentations
        aud.play('short_low');
        for itic = 1 : length(timecourse)
            if trial_valid
                if pos == 1 % first is target
                    win.draw_additive(texGrat0, michelson_to_scale(c,bg_color) * timecourse(itic),[pos0,pos0+size1px]);
                    win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                else
                    win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                end
                
                win.flip();
            end
        end
        
        % ISI
        for itic = 1 : nframes_iti
            if trial_valid
                win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                win.flip();
            end
        end
        aud.play('short_high');
        for itic = 1 : length(timecourse)
            if trial_valid
                if pos == 2 % second is target
                    win.draw_additive(texGrat0, michelson_to_scale(c,bg_color) * timecourse(itic),[pos0,pos0+size1px]);
                    win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                else
                    win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
                end
                win.flip();
            end
        end
        
        % Response interval
        list_wait.start();
        for itic = 1 : (n_frames_resp-1)
            win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
            win.flip();
        end
        response = list_wait.stop();
        if pos ==1
            aud.play('long_low');
        else
            aud.play('long_high');
        end
        win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
        win.flip; % to play the audio, frame stolen from response interval
    end

%% pause (show only fixation cross)
    function pause_trial
        list_stop.start();
        win.draw(cross_tex0,1,[pos0-sizeCross,pos0+size0cross]);
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
%% Fixation cross drawing (coppied from Sch??tz et al.)
    function draw_Fixation
        Screen('FillOval', win.h, colorOval, [cx-d1/2 * pxPerDeg(1), cy-d1/2 * pxPerDeg(1), cx+d1/2 * pxPerDeg(1), cy+d1/2 * pxPerDeg(1)], d1 * pxPerDeg(1));
        Screen('DrawLine', win.h, colorCross, cx-d1/2 * pxPerDeg(1), cy, cx+d1/2 * pxPerDeg(1), cy, d2 * pxPerDeg(1));
        Screen('DrawLine', win.h, colorCross, cx, cy-d1/2 * pxPerDeg(1), cx, cy+d1/2 * pxPerDeg(1), d2 * pxPerDeg(1));
        Screen('FillOval', win.h, colorOval, [cx-d2/2 * pxPerDeg(1), cy-d2/2 * pxPerDeg(1), cx+d2/2 * pxPerDeg(1), cy+d2/2 * pxPerDeg(1)], d2 * pxPerDeg(1));
    end


%% save
    function save_data
        if ~exist(datadir,'dir')
            mkdir(datadir);
        end
        if exist('ProbandName','var')
            fname = [ProbandName '_' num2str(typeTemporal) '_' sprintf('%.2f',f) '_' fileTimestamp '.mat'];
        else
            fname = ['test' '_' num2str(typeTemporal) '_' sprintf('%.2f',f) '_' fileTimestamp '.mat'];
        end
        if exist('results','var')
            save(fullfile(datadir, fname), 'sampling_struct','results');
        else
            save(fullfile(datadir, fname), 'sampling_struct');
        end
        fprintf('Data saved in: %s\n', fname);
    end
%% wrap up

    function wrap_up
        save_data
        try
            win.clear;
            win.delete;
            clear win
        catch
            Screen('CloseAll')
        end
    end
end
