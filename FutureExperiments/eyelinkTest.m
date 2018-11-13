% Eyelink tests

ProbandName = 'Test';
datadir       = '/home/data/heiko/parasol';
%datadir       = '/home/heiko/MATLAB/test';
fileTimestamp = standard_now;
eye_data_path = fullfile(datadir, '/eye_data_files/');
session = 1;
eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',ProbandName,'_session_', num2str(session), '.edf'));
while exist(eyedata_fname,'file')
    session = session+1;
    eyedata_fname = fullfile(eye_data_path,strcat('experiment-1_sub_',subj,'_session_', num2str(session), '.edf'));
end
% the remote file is stored on MS-DOS, so has filename length restrictions:
edf_file_remote = strcat('s_',ProbandName(1:3),'_', num2str(session), '.edf');
%
try
    calib = load('/home/data/calibration/lcd_gray/lcd_gray2018_03_13_1620.mat');
    clut = spline(calib.measurement, calib.input, linspace(0,1,(2^14))');
    clut(1:10)  = 0;
    clut(clut<0)= 0;
    
    win = window('lcd_gray', 'bg_color', .5, 'clut', clut);
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
    
    
    Eyelink('Openfile', edf_file_remote);
    
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    pause(.1)
    Eyelink('StartRecording');
    Eyelink('Message', ['start_trial']);
    pause(1)
    Eyelink('Message', ['end_trial']);
    Eyelink('StopRecording');
    
    
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
    catch e
        fprintf('Problem receiving data file ''%s''\n', edf_file_remote );
        error(e)
    end
    % Shutdown Eyelink:
    Eyelink('Shutdown');
    
catch eyelink_error
    % Shutdown Eyelink:
    Eyelink('Shutdown');
    % Close window:
    sca;
    % Restore keyboard output to Matlab:
    ListenChar(0);
    commandwindow;
    disp(eyelink_error);
end