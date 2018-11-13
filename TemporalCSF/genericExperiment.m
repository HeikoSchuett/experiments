function genericExperiment(frequencies, phases, contrasts, ...
    NperContrast,ProbandName,results)
%function mask_trials(frequencies, phases, thetas, contrasts, NperContrast)
% here NperContrast is the number of trials for each of the frequencies
% frequencies(1) is the primary frequency
% contrasts should be a length(frequencies)xn array to give contrast values
% for each frequency seperately
%
% in this version the gratings are presented very shortly in pseudo random
% order and with a 1D 1/f noise mask.
%% Setup parameters

rng('shuffle');

if nargin < 5
    NperContrast(1)=90;
    if length(frequencies)>1
        for iFreq=2:length(frequencies)
            NperContrast(iFreq) = 5;
        end
    end
end


img_size = [512 512];
bg_color = 0.5;

% How many trials per feedback screen
n_feedback=20;

aud_volume = 0.5;
%contrasts = michelson_to_scale(contrasts, bg_color);

if any(bg_color + contrasts > 1)
    error('too much')
end
% Trial durations in ms
%ms_stim = 5/72*1000; %now in frames!
ms_resp = 800;
%ms_isi = 250; %now in frames!
ms_iti = 700;

datadir       = '/home/data/heiko/short';
%datadir       = '/home/heiko/MATLAB/test';
fileTimestamp = standard_now;



dist_monitor = 1650; % mm
%degrees of visual angle of the stimulus
%CRT
deg_pic = 2*atan(.5*img_size(1)*.25/dist_monitor)/pi*180;
%LCD
%deg_pic=2*atan(.5*img_size(1)*.252/dist_monitor)/pi*180;



%% Dirty gamma-correction

% this file should contain the calibration data
calib = load('/home/data/calibration/crt_gray/crt_gray2013_09_26_1641.mat');

clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
clut(1:10)  = 0;
clut(clut<0)= 0;

%% start the hardware

% CRT-initialization
%win = window('crt_gray', 'bg_color', bg_color, 'clut', clut);
%aud = dpixx_audio_port('volume', aud_volume);
%list_wait = listener_buttonbox;
%list_stop = listener_buttonbox('does_interrupt', true);

%
% %Datapixx-initialization:
% win = window('viewpixx_gray', 'bg_color', bg_color, 'clut', clut);
% aud = dpixx_audio_port('volume', aud_volume);
% list_wait = listener_buttonbox;
% list_stop = listener_buttonbox('does_interrupt', true);

%Debugging at office initialization
win = window('debug', 'rect', [0 0 1150 1150], 'bg_color', bg_color);
aud = local_audio_port('volume', aud_volume);
list_wait = listener_keyboard;
list_stop = listener_keyboard('does_interrupt', true);


%% Calculate framecount per trial

%n_frames_stim = ms_stim * win.framerate / 1000;
n_frames_stim = 5;
%n_frames_mask = ms_mask * win.framerate / 1000;
n_frames_resp = ms_resp * win.framerate / 1000;
%n_frames_isi = ms_isi * win.framerate / 1000;
n_frames_isi = 44;
n_frames_iti = ms_iti * win.framerate / 1000;


%% Define the timecourse
timecourse = ones([n_frames_stim,1]);

%% Create textures

textures = cell(length(frequencies), length(phases),2);
for i = 1 : length(frequencies)
    f = frequencies(i)*deg_pic;
    for j = 1 : length(phases)
        p = phases(j);
        %vertical = 1
        grating = grating_sine(...
            'size', img_size, ...
            'frequency', f, ...
            'phase', p, ...
            'theta', 0);
        textures{i, j, 1} = win.make_texture(grating);
        % horizontal = 2
        grating = grating_sine(...
            'size', img_size, ...
            'frequency', f, ...
            'phase', p, ...
            'theta', 0.25);
        textures{i, j, 2} = win.make_texture(grating);
    end
end
cross_grate = patch_crosshair('size', img_size, ...
    'length', 40, 'color', .4,'thickness',5);
cross_tex = win.make_texture(cross_grate.patch);
% win_grate = distrib_tapered_cosine('size', img_size, ...
%                                    'alpha', .3);
window_tex = win.make_mask_texture('tapered_cosine', ...
    'size', img_size, ...
    'alpha', .5);





%% Generate sounds

aud.create_beep('short_low' , 'low' , 5/72, 0.25);
aud.create_beep('short_med' , 'med' , 5/72, 0.25);
aud.create_beep('short_high', 'high', .1, 0.25);
aud.create_beep('long_low'  , 'low' , .4, 0.25);
aud.create_beep('long_high' , 'high', .4, 0.25);


%% Define a trial

    function response = trial(f, c, pos)
        %function trial(f, c, pos)
        %
        % Run one trial
        
        % ITI
        for itic = 1 : n_frames_iti
            win.draw(cross_tex);
            win.flip();
        end
        
        % First interval
        aud.play('short_low');
        for itic = 1 : length(timecourse)
            if pos == 1
                win.draw_additive(f, c * timecourse(itic));
                win.draw_mask(window_tex);
            end
            win.draw(cross_tex);
            win.flip();
        end
        
        
        
        
        % ISI
        for itic = 1 : n_frames_isi
            win.draw(cross_tex);
            win.flip();
        end
        
        % Second interval
        aud.play('short_high');
        for itic = 1 : length(timecourse)
            if pos == 2
                win.draw_additive(f, c * timecourse(itic));
                win.draw_mask(window_tex);
            end
            win.draw(cross_tex);
            win.flip();
        end
        
        
        
        % Response interval
        list_wait.start();
        for itic = 1 : n_frames_resp
            win.draw(cross_tex);
            win.flip();
        end
        response = list_wait.stop();
        
    end

%% Define the clean-up function

    function wrap_up
        win.flip();
        if exist('ProbandName','var')
            fname = [ProbandName fileTimestamp];
            save(fullfile(datadir, fname), 'results', 'frequencies', 'contrasts');
        else
            fname = ['intermixedTrials' fileTimestamp];
            save(fullfile(datadir, fname), 'results', 'frequencies', 'contrasts');
        end
        Screen('Close')
        disp(sprintf('Data saved in: %s', fname));
        
        plot_experiment(results);
        disp('   contrast, Ncorrect,    N');
        disp(accumulate_trials(results(:,2:4)));
        
    end

%% Define pause trial with future target shown
    function pause_grating(stim_tex, scale)
        list_stop.start();
        win.draw_additive(stim_tex, scale);
        win.draw_mask(window_tex);
        win.draw(cross_tex);
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
%% pseudo randomization-> always a primary frequency grating after 2 consecutive trials of secondaries

    function results=pseudo_randomization(results)
        max_non_primary=2; % How many trials of non primary in a row allowed >=1 !
        
        results=results(randperm(size(results,1)),:);
        if length(frequencies) > 1
            not_accepted=true;
            primary=frequencies(1);
            while not_accepted
                not_accepted=false;
                for iline=1:(size(results,1)-max_non_primary)
                    if all(results(iline:(iline+max_non_primary),1)~=primary)
                        who_jumps=randi(max_non_primary+1)-1;
                        where_to_jump=randi(size(results,1));
                        % change rows
                        line=results(iline+who_jumps,:);
                        results(iline+who_jumps,:)=results(where_to_jump,:);
                        results(where_to_jump,:)=line;
                        not_accepted=true;
                    end
                end
            end
        end
    end
%%  true randomization
    function results = randomization(results)
        results=results(randperm(size(results,1)),:);        
    end

%% Run the experiment
%prepare list of trials to do:
%results=[frequency,contrast,truePos,response,phase,rtime]

if ~exist('results','var')
    freq_list=[];
    contrast_list=[];
    for iFreq=1:length(frequencies)
        freq_list=[freq_list;repmat(frequencies(iFreq),[NperContrast(iFreq)*size(contrasts,2),1])];
        for icontrast=1:size(contrasts,2)
            contrast_list=[contrast_list;repmat(contrasts(iFreq,icontrast),[NperContrast(iFreq),1])];
        end
    end
    results=[freq_list,contrast_list,NaN(length(freq_list),4)];
    %results=results(randperm(size(results,1)),:);
    results=pseudo_randomization(results);
    firstTrial=1;
else
    warning('You continue an older experiment. We did not test anything about it')
    firstTrial=find(isnan(results(:,4)),1);
    if isempty(firstTrial)
        error('Experiment already finished')
    end
end
% results = zeros(NperContrast*length(contrasts)*length(frequencies),6);
% zw=repmat(frequencies,[length(contrasts)*NperContrast,1]);
% results(:,1)=zw(:);
% zw=repmat(contrasts,[NperContrast,1]);
% results(:,2)=repmat(zw(:),[length(frequencies),1]);
try
    % Wait untill subject is ready
    stim_tex = textures{1, randi(length(phases)),2};
    pause_grating(stim_tex,3*michelson_to_scale(max(contrasts(:)),bg_color));
    % Loop through textures and contrasts
    for itrial=firstTrial:size(results,1)
        % Run the trial
        pos = randi(2);
        stim_phase_pos = randi(length(phases));
        results(itrial,5)=phases(stim_phase_pos);
        stim_tex = textures{find(frequencies==results(itrial,1),1), stim_phase_pos,2}; % 2 for horizontal
        results(itrial,3) = pos;
        resp = trial(stim_tex, michelson_to_scale(results(itrial,2),bg_color), pos);
        
        % Save the responses
        [press, rt] = resp.get_presses('first');
        results(itrial,4) = press;
        results(itrial,6) = rt;
        
        
        % every n_feedback trial take a brake & save & give feedback
        if mod(itrial,n_feedback)==0 
            if exist('ProbandName','var')
                fname = [ProbandName fileTimestamp];
                save(fullfile(datadir, fname), 'results', 'frequencies', 'contrasts');
            else
                fname = ['intermixedTrials' fileTimestamp];
                save(fullfile(datadir, fname), 'results', 'frequencies', 'contrasts');
            end
            % calculate performance
            correct=results(:,3)==results(:,4);
            n_correct_total=sum(correct);
            n_total=itrial;
            n_correct_last=sum(correct((itrial-n_feedback+1):itrial));
            % give feedback
            win.pause_trial(list_stop,[num2str(itrial), ' trials passed out of ', num2str(size(results,1)),...
                '\n\n overall ', num2str(100 *n_correct_total/n_total,'%2.1f'), ' % correct'...
                '\n\n ', num2str(n_correct_last), ' of the last ', num2str(n_feedback), ' Trials correct']);
            % Wait untill subject is ready
            if ~mod(itrial,360)==0 && itrial~=size(results,1)
                stim_tex = textures{1, randi(length(phases)),2};
                pause_grating(stim_tex,3*michelson_to_scale(max(contrasts(:)),bg_color));
            end
        end
        
        % forced break
        if mod(itrial,360)==0 && itrial~=size(results,1)
            for i=1:5
                win.draw_text(sprintf('you should take a short break now! \n %d of 5 minutes passed',(i-1)));
                win.flip;
                pause(60);
            end
            stim_tex = textures{1, randi(length(phases)),2};
            pause_grating(stim_tex,3*michelson_to_scale(max(contrasts(:)),bg_color));
        end
        if itrial==size(results,1)
            win.pause_trial(list_stop, 'You made it! \n Press any button to end');
        end
    end
catch e
    wrap_up();
    rethrow(e);
end

wrap_up();

end