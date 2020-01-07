module constants_mod
  implicit none
  integer, parameter :: DP=kind(1.0d0)
  real(DP), parameter :: PI=3.1415926535897932384626433832795_DP
  complex(DP), parameter :: zi=(0.0_DP, 1.0_DP)
  complex(DP), parameter :: sig=(1.0_DP, 0.0_DP)
  complex(DP), parameter :: lam=(0.333333333333333_DP, 0.0_DP)
  complex(DP), parameter :: exh=(1.0_DP, 1.0_DP)
  public
end module
