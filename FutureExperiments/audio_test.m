aud.h = PsychPortAudio('Open');
PsychPortAudio('Volume',aud.h,0.5);
aud.sample_rate = 9800;

aud_short_low = create_beep(220, 5/72, 0.25);
aud_short_high = create_beep(1000, .1, 0.25);

PsychPortAudio('FillBuffer', aud.h, [aud_short_low; aud_short_low]);
PsychPortAudio('Start', aud.h, 1, 0 , 1);

WaitSecs(1);

PsychPortAudio('FillBuffer', aud.h, [aud_short_high; aud_short_high]);
PsychPortAudio('Start', aud.h, 1, 0 , 1);

PsychPortAudio('Close');