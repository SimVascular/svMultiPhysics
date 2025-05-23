c Value of your unknown at time equal to zero (This is going to be used
c ONLY when you are initiating simulation) 

c     Variable definitions
c     x(1) = p_LV       (left ventricle pressure)
c     x(2) = p_LV_cap   (left ventricle pressure on cap (not used, just need it to compute resistance for cap))
c     x(3) = V_RV       (right ventricle volume)
c     x(4) = V_LA       (left atrium volume)
c     x(5) = V_RA       (right atrium volume)
c     x(6) = p_AR_SYS   (arterial systemic pressure)
c     x(7) = p_VEN_SYS  (venous systemic pressure)
c     x(8) = p_AR_PUL   (arterial pulmonic pressure)
c     x(9)= p_VEN_PUL  (venous pulmonic pressure)
c     x(10)= Q_AR_SYS   (arterial systemic flow rate)
c     x(11)= Q_VEN_SYS  (venous systemic flow rate)
c     x(12)= Q_AR_PUL   (arterial pulmonic flow rate)
c     x(13)= Q_VEN_PUL  (venous pulmonic flow rate)
c
c     Auxillary variables
c     x(14)= V_LV       (left ventricle volume)
c     
c     Note:
c     Q(1) = LV flow rate (out of LV)
c     Q(2) = LV cap flow rate


      tZeroX = (/ 
     &  5.6512675762742335,    ! p_LV_0      (mmHg)
     &  5.6512675762742335,    ! p_LV_cap_0    (mmHg)
     &  128.58981029386334,      ! V_RV_0        (mL)
     &  58.761096293979925,    ! V_LA_0        (mL)
     &  76.8340776729488,    ! V_RA_0        (mL)
     &  63.469005164828886,    ! p_AR_SYS_0    (mmHg)
     &  23.489460083645984,    ! p_VEN_SYS_0   (mmHg)
     &  15.39766270405246,    ! p_AR_PUL_0    (mmHg)
     &  12.990389112964845,    ! p_VEN_PUL_0   (mmHg)
     &  91.00177508885831,    ! Q_AR_SYS_0    (mL/s)
     &  112.86832799795421,    ! Q_VEN_SYS_0   (mL/s)
     &  75.19549067009953,    ! Q_AR_PUL_0    (mL/s)
     &  196.2167628991455,     ! Q_VEN_PUL_0   (mL/s)
     &  64.66340579470425      ! V_LV_0        (mL)
     & /)
    




