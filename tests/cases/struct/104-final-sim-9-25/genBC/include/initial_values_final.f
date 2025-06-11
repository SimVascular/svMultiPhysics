c Value of your unknown at time equal to zero (This is going to be used
c ONLY when you are initiating simulation) 
c
c     x(1) = p_v      (ventricle pressure)
c     x(2) = p_v_cap  (ventricle pressure on cap (not used, just need it to compute resistance for cap))
c     x(3) = p_p      (proximal pressure)
c     x(4) = q_p      (proximal flowrate)
c     x(5) = p_d      (distal pressure)
c     x(6) = V_v      (ventricle volume)
c
c Pressures in Pa, flowrates in m^3/s. Values from Martin's paper
      tZeroX = (/ 8.0*(133.322),
     &   8.0*(133.322), 
     &   61.8*(133.322),
     &   38.3*(1e-06),
     &   59.7*(133.322),
     &   0.0 /)
c      tZeroX = (/ 0.0*(133.322),
c     &   0.0*(133.322), 
c     &   0.0*(133.322),
c     &   0.0*(1e-06),
c     &   0.0*(133.322),
c     &   0.0 /)

