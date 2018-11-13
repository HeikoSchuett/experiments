function spectra = correct_spectra(spectra)
% function spectra = correct_spectra(spectra)
% as the spectrometer so far returns the old spectrum when it was too dark
% we correct it here to the lowest possible spectrum
%
% If we someday get a measurement of the low light levels, we might use it
% to do this correction (e.g. if we get the Konica minolta)

for iRGB = 1:4
    increase = true;
    i = 0;
    while increase
        i = i+1;
        if ~all(spectra{iRGB}(:,i+1)==spectra{iRGB}(:,i+2))
            increase = false;
        end
    end
    if i > 1
        spectra{iRGB}(:,2:(i+1)) = repmat(spectra{iRGB}(:,i+2),1,i);
    end
end