clear ; 
T             = 2001 ;                                     % number of time bins in sample
t             = 0:1:(T-1) ;                                % time vector in units of bins
acorr         = exp(-t/50).*cos(t/20) ;                    % hypothetical autocorrelation function
                                                           % replace this with the actual autocorrelation of the neuron (one sided) 
C             = toeplitz(acorr) ;                          % converting the autocorrelation function to a covariance matrix 
x             = mvnrnd(zeros(1,T),C) ;                     % sampling a Gaussian with a specified autocorrelation function
[~,i_x]       = sort(x,'ascend') ;                         % finding the rank of each value of x 
r_emp         = random('Exponential',1,1,1000) ;           % hypothetical empirical distribution of firing rates 
                                                           % replace this with the actual firing rates of the neuron 
rx_match(i_x) = prctile(r_emp,linspace(0,100,length(x))) ; % matching the histograms of x and the empirical distribution of firing rates
                                                           % this gives a time series of firing rates that have the desired distribution
                                                           % and approximately the desired autocorrelation function.
                                                           % The autocorrelation matches only approximately because the distribution of 
                                                           % x is Gaussian, and the firing rate distribution is not Gaussian
spks_match    = poissrnd(rx_match) ;                       % sampling spikes from a Poisson process with the matched firing rates