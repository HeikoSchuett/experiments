function calibImages
% function calibImages
% this function displays some basic calibration images on the screen, which
% can be ended by pressing any key on the keyboard

try
    close all;
    clear all;
    sca
    
    PsychImaging('PrepareConfiguration');
    PsychDefaultSetup(0)
    PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows')
    % the following three should not be used in any experiment
    Screen('Preference', 'VisualDebugLevel', 3);
    Screen('Preference', 'SuppressAllWarnings', 1);
    Screen('Preference', 'SkipSyncTests', 1);
    
    screenNumber = 0;
    
    
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    
    
    if screenNumber ==0
        [windowNr,rect]=PsychImaging('OpenWindow',screenNumber,black+0.5*(white-black),[100 100 400 400]);
    else
        [windowNr,rect]=PsychImaging('OpenWindow',screenNumber,black+0.5*(white-black));
    end
    LoadIdentityClut(windowNr);
    KbStrokeWait;
    %% zeroth image: one pixel wide border of screen
    Screen('FillRect', windowNr, white, [rect(1),rect(2),rect(1)+1,rect(4)]);
    Screen('FillRect', windowNr, white, [rect(1),rect(2),rect(3),rect(2)+1]);
    Screen('FillRect', windowNr, white, [rect(3)-1,rect(2),rect(3),rect(4)]);
    Screen('FillRect', windowNr, white, [rect(1),rect(4)-1,rect(3),rect(4)]);
    Screen('Flip',windowNr)
    KbStrokeWait;
    %% first image: 4 rows of rectangles down the screen, gray and the 3 guns
    Nsteps = 25;
    steps = round(linspace(black,white,Nsteps));
    y     = round(linspace(rect(2),rect(4),Nsteps+1));
    xdiff = rect(3)-rect(1);
    for iCol = 1:4
        for i = 1:Nsteps
            x1 = round(0.04*xdiff+rect(1)+(iCol-1)*0.24*xdiff);
            x2 = round(0.24*xdiff+rect(1)+(iCol-1)*0.24*xdiff);
            ColorRect = [x1,y(i),x2,y(i+1)];
            switch iCol
                case 1
                    Screen('FillRect', windowNr, steps(i), ColorRect);
                case 2
                    Screen('FillRect', windowNr, steps(i)*[1,0,0], ColorRect);
                case 3
                    Screen('FillRect', windowNr, steps(i)*[0,1,0], ColorRect);
                case 4
                    Screen('FillRect', windowNr, steps(i)*[0,0,1], ColorRect);
            end
        end
    end
    Screen('Flip',windowNr)
    % Clear the screen.
    KbStrokeWait;
    sca;
catch e
    sca;
    rethrow(e)
end
