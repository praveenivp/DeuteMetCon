function faxis=calcFreqAxis(dt,samples)

faxis=(floor(-samples/2):floor(samples/2)-1)'/(dt*samples);
end