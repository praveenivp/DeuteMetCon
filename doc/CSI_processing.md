# CSI processing steps

## Phase correction
As there is a delay between the excitation pulse and ADC, we need to correct for this phase if we want to use real part of the spectrum. There are two options 
1. Use Burg's method to extrapolate till the beginning of time.
2. Add linear phase shift with respect to the frequency.
