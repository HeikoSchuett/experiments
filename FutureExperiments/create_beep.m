%% create beeps
    function wave = create_beep(frequency, duration, volume)
        n = 9800 * duration;
        f = 2 * pi * frequency;
        wave = sin(f * (1 : n) / 9800);
        wave = wave * volume;  
    end