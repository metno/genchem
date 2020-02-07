module ChemFields_mod
use ChemDims_mod,   only: NSPEC_TOT ! => No. species 
implicit none
private

  real, save, public :: cell_tinv  ! 1/T where T is in Kelvin
  real, save, public, dimension(NSPEC_TOT):: &
     x, xold ,xnew  ! Work arrays [molecules/cm3]

end module ChemFields_mod


