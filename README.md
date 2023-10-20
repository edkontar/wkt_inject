# wkt_inject
WKT code (https://github.com/edkontar/wkt) with added injection of electrons (the code starts without an electron beam)
The equations are modified to include a source of electrons 
$$\frac{\partial f}{\partial t}=.... +g(v)\exp \left(-(t-t_0)^2/\tau^2\right)/\tau $$
where $\tau$ is the characteristic injection time, $t_0$ is the maximum injection time, 
and the velocity distribution of electrons is a beam $g(v)\propto \exp(-(v-v_b)^2/v_b^2)$ centered at beam velocity $v_b$.
The injection is somewhat similar to Equations 1,2 in https://ui.adsabs.harvard.edu/abs/2008AnGeo..26.2435S
