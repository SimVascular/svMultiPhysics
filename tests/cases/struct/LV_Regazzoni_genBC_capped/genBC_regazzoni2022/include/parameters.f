c     Specifing the constants
      REAL(KIND=8), PARAMETER :: pi   = 3.14159265

c    Parameters of Regazzoni 2022 model

c     ECG Parameters
      REAL(KIND=8), PARAMETER :: BPM = 87.0                 ! beats per minute
      REAL(KIND=8), PARAMETER :: P02p = 0.040               ! seconds
      REAL(KIND=8), PARAMETER :: PR_interval = 0.182        ! seconds
      REAL(KIND=8), PARAMETER :: QRS_duration = 0.088       ! seconds
      REAL(KIND=8), PARAMETER :: QT_interval = 0.384        ! seconds

c     Electromechanical delay
      REAL(KIND=8), PARAMETER :: EMD = 0.025                ! seconds

c     Contraction systole ratio
      REAL(KIND=8), PARAMETER :: CSR = 1./4.

c     Intermediate timing parameters
      REAL(KIND=8), PARAMETER :: t_start_A = EMD                        ! seconds
      REAL(KIND=8), PARAMETER :: t_end_A = PR_interval + QRS_duration   ! seconds
      REAL(KIND=8), PARAMETER :: t_start_V = PR_interval + EMD          ! seconds
      REAL(KIND=8), PARAMETER :: t_end_V = PR_interval + QT_interval    ! seconds

      REAL(KIND=8), PARAMETER :: T_sys_A = t_end_A - t_start_A         ! seconds
      REAL(KIND=8), PARAMETER :: T_sys_V = t_end_V - t_start_V         ! seconds


c    Timing Parameters
      REAL(KIND=8), PARAMETER :: T_HB = 60.0 / BPM

c     For double cosine activation model
      REAL(KIND=8), PARAMETER :: time_C_LA = t_start_A
      REAL(KIND=8), PARAMETER :: T_C_LA = CSR * T_sys_A
      REAL(KIND=8), PARAMETER :: T_R_LA = ( 1.0 - CSR ) * T_sys_A

      REAL(KIND=8), PARAMETER :: time_C_RA = t_start_A
      REAL(KIND=8), PARAMETER :: T_C_RA = CSR * T_sys_A
      REAL(KIND=8), PARAMETER :: T_R_RA = ( 1.0 - CSR ) * T_sys_A

      REAL(KIND=8), PARAMETER :: time_C_LV = t_start_V
      REAL(KIND=8), PARAMETER :: T_C_LV = CSR * T_sys_V
      REAL(KIND=8), PARAMETER :: T_R_LV = ( 1.0 - CSR ) * T_sys_V

      REAL(KIND=8), PARAMETER :: time_C_RV = t_start_V
      REAL(KIND=8), PARAMETER :: T_C_RV = CSR * T_sys_V
      REAL(KIND=8), PARAMETER :: T_R_RV = ( 1.0 - CSR ) * T_sys_V

c     For double Hill activation model
      REAL(KIND=8), PARAMETER :: m1_LA = 1.32
      REAL(KIND=8), PARAMETER :: m2_LA = 13.1
      REAL(KIND=8), PARAMETER :: tau1_LA = 0.110 * T_HB
      REAL(KIND=8), PARAMETER :: tau2_LA = 0.180 * T_HB

      REAL(KIND=8), PARAMETER :: m1_RA = 1.32
      REAL(KIND=8), PARAMETER :: m2_RA = 13.1
      REAL(KIND=8), PARAMETER :: tau1_RA = 0.110 * T_HB
      REAL(KIND=8), PARAMETER :: tau2_RA = 0.180 * T_HB

      REAL(KIND=8), PARAMETER :: m1_LV = 1.32
      REAL(KIND=8), PARAMETER :: m2_LV = 27.4
      REAL(KIND=8), PARAMETER :: tau1_LV = 0.269 * T_HB
      REAL(KIND=8), PARAMETER :: tau2_LV = 0.452 * T_HB

      REAL(KIND=8), PARAMETER :: m1_RV = 1.32
      REAL(KIND=8), PARAMETER :: m2_RV = 27.4
      REAL(KIND=8), PARAMETER :: tau1_RV = 0.269 * T_HB
      REAL(KIND=8), PARAMETER :: tau2_RV = 0.452 * T_HB

c    Ventricular elastances
      REAL(KIND=8), PARAMETER :: E_LV_act = 13.381470566581264         ! mmHg/ml
      REAL(KIND=8), PARAMETER :: E_LV_pas = 0.09214479413505382
      REAL(KIND=8), PARAMETER :: E_RV_act = 1.2468235664040908
      REAL(KIND=8), PARAMETER :: E_RV_pas = 0.030649401619707504

c    Atrial elastances
      REAL(KIND=8), PARAMETER :: E_LA_act = 0.15
      REAL(KIND=8), PARAMETER :: E_LA_pas = 0.19549185667226282
      REAL(KIND=8), PARAMETER :: E_RA_act = 0.15
      REAL(KIND=8), PARAMETER :: E_RA_pas = 0.06736966026440389

c    Atrial rest volumes
      REAL(KIND=8), PARAMETER :: V0_LA = 26.23534455334443
      REAL(KIND=8), PARAMETER :: V0_RA = 41.680015938842274

c    Ventricular rest volumes
      REAL(KIND=8), PARAMETER :: V0_LV = 32.41857776349184
      REAL(KIND=8), PARAMETER :: V0_RV = 72.05452710344869

c    Systemic circulation
      REAL(KIND=8), PARAMETER :: R_AR_SYS = 0.4477379610775147
      REAL(KIND=8), PARAMETER :: C_AR_SYS = 1.5171714595050794
      REAL(KIND=8), PARAMETER :: L_AR_SYS = 0.005
      REAL(KIND=8), PARAMETER :: Z_AR_SYS = 0.0
      REAL(KIND=8), PARAMETER :: R_VEN_SYS = 0.18711232422995847
      REAL(KIND=8), PARAMETER :: C_VEN_SYS = 60.0
      REAL(KIND=8), PARAMETER :: L_VEN_SYS = 0.0005

c    Pulmonary circulation
      REAL(KIND=8), PARAMETER :: R_AR_PUL = 0.032
      REAL(KIND=8), PARAMETER :: C_AR_PUL = 10.0
      REAL(KIND=8), PARAMETER :: L_AR_PUL = 0.0005
      REAL(KIND=8), PARAMETER :: Z_AR_PUL = 0.0
      REAL(KIND=8), PARAMETER :: R_VEN_PUL = 0.035
      REAL(KIND=8), PARAMETER :: C_VEN_PUL = 16.0
      REAL(KIND=8), PARAMETER :: L_VEN_PUL = 0.0005

c    Valve resistances
      REAL(KIND=8), PARAMETER :: R_min = 0.005
      REAL(KIND=8), PARAMETER :: R_max = 5.0


