function the_am_noise = amplitude_modulated_noise_example()
%the_am_noise = amplitude_modulated_noise_example
% Sample script showing how to amplitude-modulate some sound.
% Along the way, read/write a wave file and phase-shuffle to make noise with
% matching spectral profile.
%
% Be careful with the amount of modulation! Don't just do whatever, you
% need to calculate how much modulation of the signal is needed for given
% modulation at the sensory and perceptual levels!
%
% DGD

% Matlab's own audio file example, loads a waveform and sampling rate
load handel

% This is how you save some data to disk as a wave file
audiowrite('handel.wav',y,Fs);

% This is how you read it
[y,Fs] = audioread('handel.wav');

% Amplitude Envelope.
switch 1
    case 1 % AE of original input
        % Before we amplitude modulate this sample let's first get its original AE
        % Some parameters
        len = .025; % Length of running window, in seconds
        window_span = round(Fs*len); % Length of running window, in samples
        
        bins = 0:window_span:length(y);
        envelope_x = zeros(length(bins)-1,1);
        for c = 2:numel(bins)
            envelope_x(c-1) = std(y(bins(c-1)+1:bins(c)));
        end
        envelope_x = interp(envelope_x,floor(length(y)/length(envelope_x))); % Upsample to the original signal
        
        % They might end up being not exactly the same length
        if size(envelope_x,1)<size(y,1)
            envelope_x = [envelope_x;zeros(length(y)-length(envelope_x),1)];
        end
        if size(envelope_x,1)>size(y,1)
            envelope_x = envelope_x(1:size(y,1));
        end
        % Here the env is not normalized
        
    case 2 % A sine wave AE
        cycles_per_second = 30; % Hz
        modulation_amplitude = .5; % amplitude is half-amplitude!
        envelope_x = sin((1:size(y,1))'*2*pi*cycles_per_second)*modulation_amplitude-modulation_amplitude+1;
        % Safer to AM down not up to avoid peaking, so norm the env to 1
end

% Noise. 1/f^beta betaE{1,2} noise is the safest. Many methods available
switch 1
    case 1 % Noise is the phase-shuffled version of the original sound
        noise = phase_randomize(y);
    case 2 % Pink noise
        noise = spatialPattern(size(y),-1);
end
noise = noise./max(abs(noise)); % Normalize it

% Modulate the noise
the_am_noise = noise.*envelope_x; % This is the core of the amplitude modulation
the_am_noise = the_am_noise.*(.9*max(abs(the_am_noise))); % Normalize to near-unit variance
the_am_noise = the_am_noise./std(the_am_noise)*std(y); % Normalize to the variance of x

% the_am_y = y.*envelope_x;
% the_am_y = the_am_y.*(.9*max(abs(the_am_y)));
% the_am_y = the_am_y./std(the_am_y)*std(y);

% Verify visually
subplot(3,1,1)
plot(y)
subplot(3,1,2)
plot(envelope_x)
subplot(3,1,3)
plot(the_am_noise)

% audiowrite('handel_shuffled.wav',y,Fs);
% sound(the_am_noise)

end


function Y = phase_randomize(X)
% Returns a phase-randomized version of the input data X.
if size(X,2)>size(X,1);X=X';end
for c=1:size(X,2)
    x=X(:,c);
    pad_flag=0;if mod(numel(x),2)==1;x=x(1:end-1);pad_flag=1;end
    n = numel(x);
    y = fft(x); % Get spectrum
    % Add random phase shifts (negative for conjugates), preserve DC offset
    rnd_theta = -pi + (2*pi).*rand(n/2-1,1);
    y(2:n/2) = y(2:n/2).*exp(1i*rnd_theta);
    y(n/2+2:n) = y(n/2+2:n).*exp(-1i*flipud(rnd_theta));
    % return phase-randomized data
    y = ifft(y);
    if pad_flag==1;y=vertcat(y,0);end
    Y(:,c)=y;
end
end