
### phase cycles
You can pick the `shift` parameter such that the first PC is passband(180 deg phase increment)
```
NRep=54;
Range=360;
Shift=180-(180/NRep);
phi_vec = (Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift
```