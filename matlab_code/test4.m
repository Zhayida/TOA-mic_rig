% tdoa (3 dual-rigs + 5 sounds, slightly overdetermined

if 1
    clearvars
    %%
    m = 3;
    n = 5;
    
    
    options.dim = 3;
    options.origin = 1;
    options.rsyn = 0;
    
    [data] = generate_mic_rig(m,n,options);
end