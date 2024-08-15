%%
dt=1e-6; %s
t=0:dt:0.01; %s

rfdur=1.4e-3;% s

RFpulse=zeros(size(t));

RFpulse((1:round(rfdur/dt)))=1;%+round(length(RFpulse)/2))=1;

faxis=linspace(-0.5/dt,0.5/dt,length(RFpulse));
spec=abs(fftshift(fft(RFpulse,[],2),2));

ind = find(min(abs(spec - max(spec)*0.9)) == abs(spec - max(spec)*0.9));

figure,
subplot(211),plot(t,RFpulse)
subplot(212),plot(faxis,abs(fftshift(fft(RFpulse,[],2),2)))
xlabel('freq'),xlim([-1 1]*1500)
hold on
plot(faxis(ind),spec(ind),'*')
diff(faxis(ind))


