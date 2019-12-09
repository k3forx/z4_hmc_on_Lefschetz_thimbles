module generate_thimble_mod
  use constants_mod
  implicit none
  type critical_point
    complex(DP) :: val, vec
    real(DP) :: kap
  end type

contains

! subroutine make_complex_symmetric_matrix(K, z)
!   implicit none
!   complex(DP), intent(out) :: K(:,:)
!   complex(DP), intent(in) :: z(:)

!   K(1,1) = sig + 3.0_DP*lam*z(1)**2

!   return
! end subroutine

subroutine set_init_cond(z, tvec, cp, t_init, ori)
  implicit none
  type(critical_point), intent(in) :: cp
  real(DP), intent(in) :: ori, t_init
  complex(DP), intent(inout) :: z, tvec

  z = cp%val + cp%vec*exp(cp%kap*t_init)*ori
  tvec = cp%vec*exp(cp%kap*t_init)*ori

  return
end subroutine

function action_val(z, ord) result (S)
  implicit none
  complex(DP), intent(in) :: z
  integer, intent(in) :: ord
  complex(DP) :: S

  if (ord == 0) then
    S = 0.5_DP*sig*z**2 + 0.25_DP*lam*z**4 + exh*z
  else if (ord == 1) then
    S = sig*z + lam*z**3 + exh
  else if (ord == 2) then
    S = sig + 3.0_DP*lam*z**2
  end if

  return
end function

subroutine runge_kutta_new_config(z1, z0, h)
  implicit none
  real(DP), intent(in) :: h
  complex(DP), intent(in) :: z0
  complex(DP), intent(out) :: z1
  complex(DP) :: k1, k2, k3, k4
  integer, parameter :: ord = 1

  k1 = h * conjg(action_val(z0, ord))
  k2 = h * conjg(action_val(z0 + 0.5_DP*k1, ord))
  k3 = h * conjg(action_val(z0 + 0.5_DP*k2, ord))
  k4 = h * conjg(action_val(z0 + k3, ord))

  z1 = z0 + (k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine

subroutine runge_kutta_new_vector(tvec1, z0, tvec0, h)
  implicit none
  real(DP), intent(in) :: h
  complex(DP), intent(in) :: z0, tvec0
  complex(DP), intent(out) :: tvec1
  complex(DP) :: k1, k2, k3, k4
  integer, parameter :: ord = 2

  k1 = h * conjg(tvec0 * action_val(z0, ord))
  k2 = h * conjg(tvec0 * action_val(z0 + 0.5_DP*k1, ord))
  k3 = h * conjg(tvec0 * action_val(z0 + 0.5_DP*k2, ord))
  k4 = h * conjg(tvec0 * action_val(z0 + k3, ord))

  tvec1 = tvec0 + (k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine

subroutine solve_flow_eq (z0, z1, tvec0, tvec1, trange)
  implicit none
  complex(DP), intent(inout) :: z0, z1, tvec0, tvec1
  real(DP), intent(in) :: trange
  real(DP) :: t, dt, w
  integer :: itre, ii

  itre = 10000
  dt = trange/real(itre, kind=DP)

  do ii = 1, itre
    !
    ! Runge kutta method (4th order)
    !
    call runge_kutta_new_config(z1, z0, dt)
    call runge_kutta_new_vector(tvec1, z0, tvec0, dt)

    !
    ! update variable
    !
    z0 = z1
    tvec0 = tvec1

    ! w = exp(-real(action_val(z0, 0)))
    ! if (w <= 1E-20) return

  end do

  return
end subroutine
end module
!  call make_complex_symmetric_matrix(K_mat,fp)
!  call takagi_factorization(kap, vec, conjg(K_mat))
