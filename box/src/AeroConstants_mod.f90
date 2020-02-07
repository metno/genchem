module AeroConstants_mod

  ! BoxChem/ESX/EMEP need some specific calculations of aerosol
  ! surface area, and we define 7 aerosol types
  ! We end up with variables here to avoid circularity in
  ! Makefile dependencies
   
  implicit none

  private 

  !A2018 - added to allow emepctm-like aerosol reactions
  !        so we can refer to eg AERO%PM_F

   integer, parameter, public :: NSAREA_DEF = 7 ! skip ORIG=Riemer
   type, public :: aero_t
     integer :: & ! 
       SIA_F=1,PM_F=2,SS_F=3,DU_F=4,SS_C=5,DU_C=6,PM=7,NSAREA=NSAREA_DEF
   end type aero_t
   type(aero_t), public, save :: AERO = aero_t()

end module AeroConstants_mod
