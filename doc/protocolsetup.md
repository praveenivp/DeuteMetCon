
### phase cycle setup
```
NRep=54;
Range=360;
Shift=ceil(180/NRep);
phi_vec = (Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift
```