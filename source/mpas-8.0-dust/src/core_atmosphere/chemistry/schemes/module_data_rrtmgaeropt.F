
MODULE module_data_rrtmgaeropt
    ! use mpas_atmchem_vars, only:nswbands, nlwbands
      IMPLICIT NONE
!jwfeng comment nswbands and nlwbands which have been defined at Registry.xml
      INTEGER nswbands,nlwbands   ! wave bands for rrtmg radiation scheme
      !jwfeng: Use all 14 sw bands of OPAC optical parameters
      PARAMETER (nswbands = 14,nlwbands=16)

!*************************************************************
!czhao hard coding the refractive index of water and aerosols
      !water
      real,dimension(1:nswbands),save ::  refrwsw,refiwsw
      real,dimension(1:nlwbands),save ::  refrwlw,refiwlw
      !data refrwsw /1.35,1.34,1.33,1.33/
      !data refiwsw /1.524e-8,2.494e-9,1.638e-9,3.128e-6/
      data refrwsw /1.40804194921731, 1.32446157174469, 1.27672583279689, 1.30184716756451, &
                    1.31205127522499, 1.32130110609241, 1.32283620723917, 1.32683645035353, &
                    1.33096787187364, 1.33366995434286, 1.33954684893077, 1.34852085514705, &
                    1.36199998855591, 1.29078301665348/
      data refiwsw /0.014204204351483, 0.157659209404304, 0.001516344010531, &
                    0.001159062473243, 0.000235527521489, 0.000171196911521, &
                    2.43626001650e-05, 3.12758360266e-06, 3.74323598365e-08, &
                    1.63841034374e-09, 2.49434955910e-09, 1.52413800466e-08, &
                    3.35000009955e-08, 0.034382551791929/
      data refrwlw /1.532,1.524,1.420,1.274,1.161,1.142,1.232,1.266,1.296, & 
      1.321,1.342,1.315,1.330,1.339,1.350,1.408/ 
      data refiwlw / 0.336,0.360,0.426,0.403,0.321,0.115,0.0471,0.039,0.034, &
      0.0344,0.092,0.012,0.013,0.01,0.0049,0.0142/

      !dust
      real,dimension(1:nswbands),save ::  refrsw_dust,refisw_dust
      real,dimension(1:nlwbands),save ::  refrlw_dust,refilw_dust
      !data refrsw_dust /nswbands*1.530/
      !data refrsw_dust /nswbands*1.550/
      !data refisw_dust /0.024,0.0135,0.0063,0.004/
      !data refisw_dust /nswbands*0.003/
      data refrsw_dust  /1.510, 1.520, 1.523, 1.529, 1.530, 1.530, 1.530, &
                         1.530, 1.530, 1.530, 1.530, 1.530, 1.530, 1.180/
      data refisw_dust  /0.0184, 0.0284, 0.0118, 0.0082, 0.0066, &
                         0.0056, 0.0051, 0.0040, 0.0040, 0.0063, &
                         0.0135, 0.0244, 0.0300, 0.1014/

      data refrlw_dust /2.340,2.904,1.748,1.508,1.911,1.822,2.917,1.557, &
                        1.242,1.447,1.432,1.473,1.495,1.500,1.500,1.510/
      data refilw_dust /0.7,0.857,0.462,0.263,0.319,0.26,0.65,0.373,0.093, &
                        0.105,0.061,0.0245,0.011,0.008,0.0068,0.018/

      !BC
      real,dimension(1:nswbands),save ::  refrsw_bc,refisw_bc
      real,dimension(1:nlwbands),save ::  refrlw_bc,refilw_bc
      data refrsw_bc /nswbands*1.95/
      data refisw_bc /nswbands*0.79/
      data refrlw_bc /nlwbands*1.95/
      data refilw_bc /nlwbands*0.79/

      !OC
      real,dimension(1:nswbands),save ::  refrsw_oc,refisw_oc
      real,dimension(1:nlwbands),save ::  refrlw_oc,refilw_oc
      !data refrsw_oc /1.53,1.53,1.53,1.52/
      !data refisw_oc /0.00776,0.005,0.00567,0.0156/
      !data refrsw_oc /nswbands*1.45/
      !data refisw_oc /nswbands*0.0/
      data refrsw_oc /1.44300699319961, 1.41999995708466, 1.41999995708466, &
                      1.41999995708466, 1.46322364269253, 1.50999999046326, 1.50999999046326, &
                      1.51959105928615, 1.52999997138977, 1.52999997138977, 1.52999997138977, &
                      1.52999997138977, 1.52999997138977, 1.1235961228882 /
      data refisw_oc /0.0057167844689993, 0.0177692296388509, &
                      0.0106021503694548, 0.0083691408588851, 0.0162124886333433, &
                      0.0219730231735957, 0.0192865604244433, 0.0156431224426743, &
                      0.0070481515576108, 0.0056650108233176, 0.0049999998882413, &
                      0.0077604362740354, 0.0299999993294477, 0.0790744021342926/
      data refrlw_oc /1.86,1.91,1.988,1.439,1.606,1.7,1.888,2.489,1.219, &
      1.419,1.426,1.446,1.457,1.458,1.455,1.443/ 
      data refilw_oc /0.5,0.268,0.185,0.198,0.059,0.0488,0.11,0.3345,0.065, &
      0.058,0.0261,0.0142,0.013,0.01,0.005,0.0057/

      !Sea-salt
      real,dimension(1:nswbands),save ::  refrsw_seas,refisw_seas
      real,dimension(1:nlwbands),save ::  refrlw_seas,refilw_seas
      !data refrsw_seas /1.51,1.5,1.5,1.47/
      !data refisw_seas /0.866e-6,7.019e-8,1.184e-8,0.00015/
      data refrsw_seas /1.48000001907349, 1.53384612000459, 1.43698923169757, &
                        1.44815433592094, 1.45000004768372, 1.46150568669356, 1.46918128392771, &
                        1.47000002861023, 1.49000000953674, 1.50000000000000, 1.50136713125929, &
                        1.50999999046326, 1.50999999046326, 1.40804248446564/
      data refisw_seas /0.00175734291533957, 0.00746153774952505, &
                        0.00295161276680147, 0.00127685543890676, 0.00079425554125163, &
                        0.00053808019739796, 0.00037535739277884, 0.00014987360459137, &
                        2.2824878297772e-07, 1.1842440250740e-08, 7.0193698152790e-08, &
                        1.8661635377763e-06, 4.9999998736894e-06, 0.0158148805727024/
      data refrlw_seas /1.74,1.76,1.78,1.456,1.41,1.48,1.56,1.63,1.4,1.43, &
      1.56,1.45,1.485,1.486,1.48,1.48 / 
      data refilw_seas /0.1978,0.1978,0.129,0.038,0.019,0.014,0.016,0.03,0.012, &
      0.0064,0.0196,0.0029,0.0017,0.0014,0.0014,0.00176/

      !Sulfate
      real,dimension(1:nswbands),save ::  refrsw_sulf,refisw_sulf
      real,dimension(1:nlwbands),save ::  refrlw_sulf,refilw_sulf
      !data refrsw_sulf /1.468,1.442,1.43,1.422/
      !data refrsw_sulf /nswbands*1.52/
      !data refisw_sulf /3*1.0e-9,1.75e-6/
      data refrsw_sulf /1.36690911341418, 1.31457691114656, 1.3579784668612, &
                        1.38030854221016, 1.39264476329731, 1.4045056456856, 1.41218124291971, &
                        1.42163201814358, 1.42796787473467, 1.4303349523759, 1.44164061147211, &
                        1.46764245536755, 1.48399996757507, 1.1641275585545/
      data refisw_sulf /0.157930075893293, 0.0566869128264954, &
                        0.002886343867360, 0.0014907128582646, 0.00053038523314, &
                        0.000102977118696, 1.619673579599e-05, 1.7512267758e-06, &
                        2.21435654637e-08, 9.999999939225e-09, 9.9999999392e-09, &
                        9.99999993923e-09, 9.999999939225e-09, 0.55113374603149/
      data refrlw_sulf /1.89,1.91,1.93,1.586,1.678,1.758,1.855,1.597,1.15, & 
      1.26,1.42,1.35,1.379,1.385,1.385,1.367/
      data refilw_sulf /0.22,0.152,0.0846,0.2225,0.195,0.441,0.696,0.695, & 
      0.459,0.161,0.172,0.14,0.12,0.122,0.126,0.158/

!*************************************************************
      !wavelength
      real, save :: wavmin(nswbands) ! Min wavelength (um) of interval
      !data wavmin /0.25,0.35,0.55,0.998/
      data wavmin / 3.077, 2.500, 2.150, 1.942, 1.626, 1.299, 1.242, &
                    0.778, 0.625, 0.442, 0.345, 0.263, 0.200, 3.846  /
      real, save :: wavmax(nswbands) ! Max wavelength (um) of interval
      data wavmax / 3.846, 3.077, 2.500, 2.150, 1.942, 1.626, 1.299, &
                    1.242, 0.778, 0.625, 0.442, 0.345, 0.263, 12.195 /
      !data wavmax/0.35,0.45,0.65,1.000/
      real,save :: wavmidsw(nswbands)
      !data wavmidsw  / 0.30e-4, 0.40e-4, 0.60e-4 ,0.999e-04 /

      real, save :: wavenumber1_longwave(nlwbands) !Longwave limits (cm-1)
      data wavenumber1_longwave /10.,350.,500.,630.,700.,820.,980.,1080.,1180.,1390.,1480.,1800.,2080.,2250.,2390.,2600./
      real, save :: wavenumber2_longwave(nlwbands) !Longwave limits (cm-1)
      data wavenumber2_longwave /350.,500.,630.,700.,820.,980.,1080.,1180.,1390.,1480.,1800.,2080., 2250.,2390.,2600.,3250./

      !mode or size bin
      integer,parameter :: maxd_amode=3
      integer,parameter :: ntot_amode=3
      integer,parameter :: maxd_bin=8
      integer,parameter :: ntot_bin=8

      !Chebychev polynomial
      !integer,parameter :: prefr=7,prefi=10
      integer,parameter :: prefr=7,prefi=7
      integer,parameter :: ncoef=50
      real,parameter :: rmmin=0.005e-4,rmmax=50.e-4  ! cm 
      real,save:: refrtabsw(prefr,nswbands)
      real,save:: refitabsw(prefi,nswbands)
      real,save:: refrtablw(prefr,nlwbands)
      real,save:: refitablw(prefi,nlwbands)
      !coefficients for parameterizing aerosol radiative properties
      !in terms of refractive index and wet radius
      real,save:: extpsw(ncoef,prefr,prefi,nswbands) !specific extinction
      real,save:: abspsw(ncoef,prefr,prefi,nswbands) !specific absorption
      real,save:: ascatpsw(ncoef,prefr,prefi,nswbands) !specific scattering
      real,save:: asmpsw(ncoef,prefr,prefi,nswbands) !asymmetry factor
      real,save:: sbackpsw(ncoef,prefr,prefi,nswbands)
      real,save:: pmom2psw(ncoef,prefr,prefi,nswbands) 
      real,save:: pmom3psw(ncoef,prefr,prefi,nswbands)
      real,save:: pmom4psw(ncoef,prefr,prefi,nswbands)
      real,save:: pmom5psw(ncoef,prefr,prefi,nswbands)
      real,save:: pmom6psw(ncoef,prefr,prefi,nswbands)
      real,save:: pmom7psw(ncoef,prefr,prefi,nswbands)
      real,save:: extplw(ncoef,prefr,prefi,nlwbands) !specific extinction
      real,save:: absplw(ncoef,prefr,prefi,nlwbands) !specific absorption
      real,save:: ascatplw(ncoef,prefr,prefi,nlwbands) !specific scattering
      real,save:: asmplw(ncoef,prefr,prefi,nlwbands) !asymmetry factor

      !now czhao use 0.45 instead of 0.40 becaues of incorrect AOD from 0.40
      !data wavmidsw  / 0.30e-4, 0.45e-4, 0.60e-4 ,0.999e-04 /
      real,save :: wavmidlw(nlwbands)
      complex, save :: crefwsw(nswbands) ! complex refractive index fro water
      complex, save :: crefwlw(nlwbands)



END MODULE module_data_rrtmgaeropt
