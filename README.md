# wkt_inject
WKT code (https://github.com/edkontar/wkt) with added injection of electrons (starting without any beam)
The equations are modified to include a source of electrons 
$$\frac{\partial f}{\partial t}=....-\left.v\left(\frac{\partial f_0}{\partial x}\right)\right|_{x=x_S}(v, t)$$
where $f_0(v, x, t)=g(v) \exp \left(-(x_s-v t)^2 / d^2\right)$
and velocity distribution of electrons is a power-law $g(v)=A v^{-\delta}$.
