function data = downsample_following(data, ds_target)

%sf_check = input('Is the sampling frequency 44.1 kHz? Type 1 for yes.')

sf = 44100;
ds_factor = round(sf/ds_target);
data = downsample(data,ds_factor);
    
