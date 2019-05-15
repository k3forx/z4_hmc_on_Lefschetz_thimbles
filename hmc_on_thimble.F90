module sub_hmc_on_thimble
  implicit none
  integer, parameter :: DP=kind(1.0d0)
  real(DP),parameter :: PI = 3.141592653589793_DP
  complex(DP), parameter :: zi=(0.0_DP,1.0_DP)
  complex(DP), parameter :: sig=(1.0_DP, 0.0_DP)
  complex(DP), parameter :: lam=(0.333333333333333_DP, 0.0_DP)
  complex(DP), parameter :: exh=(1.0_DP, 1.0_DP)

contains

function generate_gauss_rand() result(grand)
  use mt19937
  implicit none
  real(DP) :: grand, rr, tt

  rr = grnd()
  tt = grnd()
  rr = 1.0_DP - rr
  tt = 1.0_DP - tt
  grand = sqrt(-2.0_DP*log(rr))*cos(2.0_DP*PI*tt)

  return
end function

subroutine set_init_param(t_init, ori)
  implicit none
  real(DP), intent(out) :: t_init, ori(:)
  real(DP) :: norm
  real(DP), allocatable :: grand(:)
  integer :: ii, ndim

  ndim=size(ori(:))
  allocate (grand(1:ndim))

  norm = 0.0_DP
  do ii = 1, ndim
    grand(ii) = generate_gauss_rand()
    norm = norm + grand(ii)**2
  end do

  do ii = 1, ndim
    ori(ii) = grand(ii)*sqrt(real(ndim, kind=DP)/norm)
  end do

  t_init = -2.0_DP

  deallocate (grand)
  return
end subroutine

function calc_config(t, ori, fp, vec, kap) result(z)
  implicit none
  real(DP), intent(in) :: t, ori(:), kap(:)
  complex(DP) :: fp, vec(:)
  complex(DP) :: z, sum_grad
  integer :: ii, ndim

  ndim = size(ori(:))
  sum_grad = (0.0_DP, 0.0_DP)

  do ii = 1, ndim
    sum_grad = sum_grad + vec(ii)*exp(kap(ii)*t)*ori(ii)
  end do

  z = fp + sum_grad

  return
end function

subroutine calc_tangent_vector(tvec, vec, kap, t)
  implicit none
  complex(DP), intent(out) :: tvec(:)
  complex(DP), intent(in) :: vec(:)
  real(DP), intent(in) :: kap(:), t
  integer :: ii, ndim

  ndim = size(vec(:))
  do ii = 1, ndim
    tvec(ii) = vec(ii)*exp(kap(ii)*t)
  end do

  return
end subroutine


end module

program hmc_on_thimble
  use mt19937
  use sub_hmc_on_thimble
  implicit none
  real(DP), allocatable :: ori(:)
  real(DP) :: t
  complex(DP) :: z0, z1
  complex(DP), allocatable :: tvec0(:), tvec1(:)
  real(DP) :: p0, p1

  complex(DP) :: fp
  complex(DP), allocatable :: vec(:)
  real(DP), allocatable :: kap(:)

  allocate (tvec0(1), tvec1(1), vec(1), kap(1), ori(1))

  fp = (-9.7492178167965704E-01, -5.3953762329131794E-01)
  vec(1) = (-9.7492178167965704E-01, -5.3953762329131794E-01)
  kap(1) = 1.9647512804566389E+00


  call sgrnd(1)
  call set_init_param(t, ori)

  z0 = calc_config(t, ori, fp, vec, kap)
  call calc_tangent_vector(tvec0, vec, kap, t)

  stop
end program
