function imgnoiseout = SynCamNoises(imgsize,imgmin,baseline)
qe = 0.82;
% pxsize = 6.5;
sensitivity = 5.88;
dark_noise = 0.06 * sensitivity;
% baseline = 100;

noise_poisson = poissrnd(imgmin/sensitivity/qe,imgsize) * qe;
noise_normal = normrnd(0,dark_noise,imgsize);

imgnoiseout = uint16(max(0,(noise_poisson + noise_normal ) * sensitivity + baseline));
end