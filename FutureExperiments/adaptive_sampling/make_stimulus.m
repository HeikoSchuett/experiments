function stimulus = make_stimulus(im_size,deg_size,s1,s2,background)
% function stimulus = make_stimulus(im_size,deg_size,s1,s2,background)
% This function creates images to be used in an peripheral contrast
% detection experiment. 
% s1 and s2 are the specifications of the targets in the following format: 
% [sf,contrast,orientation,phase,std,excentricity]


%stimulus = zeros(im_size);

gab1  = gabor(s1,im_size,deg_size);
gab2  = gabor(s2,im_size,deg_size);

if isnumeric(background)
    stimulus = background + gab1 + gab2;
elseif strcmp(background,'blank')
    stimulus = 0.5 + gab1 + gab2;
elseif strcmp(background,'noise')
    stimulus = noise(im_size,deg_size,.1 ) + gab1 + gab2;
end


stimulus(:,ceil(im_size(2)/2))=0;
stimulus(ceil(im_size(1)/2),ceil(im_size(2)/2)+(-10:10))=0;
end

function gab = gabor(s,im_size,deg_size)
% generates a gabor patch
x = (1:im_size(2))./im_size(2).*deg_size(2);
y = (1:im_size(1))./im_size(1).*deg_size(1);

x = x-mean(x)-s(6); % in x dimension move out to periphery from the mean
y = y-mean(y);      % in y dimension centered

[xx,yy]= meshgrid(x,y');

xx2 = cos(s(3))*xx-sin(s(3))*yy;
yy2 = cos(s(3))*yy+sin(s(3))*xx;

gab = s(2)* exp(-(xx2.^2+yy2.^2)/(s(5)^2)).*sin(2*pi*xx2*s(1)-s(4));

end

function n = noise(im_size,deg_size,stdev)
n = randn(im_size);
m = mean(n(:));
n = n-m;
n = fft2(n);
x = 1:im_size(2);  % for frequency space
y = 1:im_size(1);
x = (x-ceil(mean(x)))./deg_size(2);  % normalize to mean in middle, and freq in cyc/deg
y = (y-ceil(mean(y)))./deg_size(1);

freq = bsxfun(@(x,y) sqrt(x.^2+y.^2),x,y');
n    = n./(.001+fftshift(freq));

n    = ifft2(n);

n = 0.5+ stdev*(n-mean(n(:)))/std(n(:));
end