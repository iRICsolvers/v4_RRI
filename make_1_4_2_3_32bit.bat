del *.obj *.mod
ifort /c /Qopenmp /O3 ./1.4.2.3/RRI_Mod.f90
ifort /c /Qopenmp /O3 ./1.4.2.3/RRI_Mod2.f90
ifort /c /Qopenmp /O3 ./1.4.2.3/RRI_Mod_Dam.f90
ifort /c /Qopenmp /O3 ./1.4.2.3/RRI_Mod_Tecout.f90
ifort /c  /Qopenmp /O3  ./1.4.2.3/RRI.f90 ./1.4.2.3/RRI_Bound.f90 ./1.4.2.3/RRI_Dam.f90 ./1.4.2.3/RRI_Div.f90  ./1.4.2.3/RRI_DT_Check.f90 ./1.4.2.3/RRI_Evp.f90 ./1.4.2.3/RRI_GW.f90 ./1.4.2.3/RRI_Infilt.f90 ./1.4.2.3/RRI_Read.f90 ./1.4.2.3/RRI_Riv.f90 ./1.4.2.3/RRI_RivSlo.f90 ./1.4.2.3/RRI_section.f90 ./1.4.2.3/RRI_Slope.f90 ./1.4.2.3/RRI_Sub.f90 ./1.4.2.3/RRI_Tecout.f90 ./1.4.2.3/RRI_TSAS.f90
ifort /MT /Qopenmp /O3 *.obj /o ../bin32/0_rri_1_4_2.exe
editbin /stack:40000000 ../bin32/0_rri_1_4_2.exe
del *.obj *.mod
