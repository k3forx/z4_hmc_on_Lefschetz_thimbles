module sub_generate_thimble
  use constants_mod
  implicit none
  type critical_point
    complex(DP) :: val, vec
    real(DP) :: kap
  end type

contains

subroutine make_complex_symmetric_matrix(K, z)
  implicit none
  complex(DP),intent(inout) :: K(:,:)
  complex(DP),intent(in) :: z(:)

  K(1,1) = sig + 3.0_DP*lam*z(1)**2

  return
end subroutine

!subroutine takagi_factorization(diag, vec, K)
!  use forpy_mod
!  implicit none
!  complex(DP),intent(in) :: K(:,:)
!  real(DP),dimension(:,:),pointer,intent(out) :: diag
!  complex(DP),dimension(:,:),pointer,intent(out) :: vec

!  integer :: ierr
!  type(module_py) :: python_mod
!  type(ndarray) :: nd_mat, nd_diag, nd_vec
!  type(tuple) :: arg
!  type(list) :: paths
!  type(object) :: res_takagi, attr

!  ierr = forpy_initialize()
!  ierr = get_sys_path(paths)
!  ierr = paths%append(".")

!  ierr = import_py(python_mod, "python_mod")

!  ierr = ndarray_create(nd_mat, K)
!  ierr = tuple_create(arg,1)
!  ierr = arg%setitem(0, nd_mat)

!  !
!  ! Takagi factorization via Python module
!  !
!  ierr = call_py(res_takagi, python_mod, "Takagi", arg)

!  !
!  ! get diagonal matrix from 'res_takagi' object
!  !
!  ierr = res_takagi%getattribute(attr, "diag")
!  ierr = cast(nd_diag, attr)
!  ierr = nd_diag%get_data(diag)
!  call attr%destroy

!  !
!  ! get orthonormal vector from 'res_takagi' object
!  !
!  ierr = res_takagi%getattribute(attr, "ortho_vec")
!  ierr = cast(nd_vec, attr)
!  ierr = nd_vec%get_data(vec)
!  call attr%destroy

!  call python_mod%destroy
!  call nd_mat%destroy
!  call nd_diag%destroy
!  call nd_vec%destroy
!  call arg%destroy
!  call paths%destroy
!  call res_takagi%destroy

!  call forpy_finalize()

!  return
!end subroutine

subroutine set_init_cond(z, tvec, cp, t_init, ori)
  implicit none
  type(critical_point), intent(in) :: cp
  real(DP), intent(in) :: ori, t_init
  complex(DP), intent(inout) :: z, tvec

  z = cp%val + cp%vec*exp(cp%kap*t_init)*ori
  tvec = cp%vec*exp(cp%kap*t_init)*ori
  write(*,'("init_cond: ",10ES24.15)') z, tvec

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
  integer,parameter :: ord = 1

  ! Generate J_{\sigmal}
  k1 = conjg(action_val(z0, ord))
  k2 = conjg(action_val(z0 + 0.5_DP*h*k1, ord))
  k3 = conjg(action_val(z0 + 0.5_DP*h*k2, ord))
  k4 = conjg(action_val(z0 + h*k3, ord))

  ! Generate K_{\sigma}
  ! k1 = -1.0_DP * conjg(action_val(z0, ord))
  ! k2 = -1.0_DP * conjg(action_val(z0 + 0.5_DP*h*k1, ord))
  ! k3 = -1.0_DP * conjg(action_val(z0 + 0.5_DP*h*k2, ord))
  ! k4 = -1.0_DP * conjg(action_val(z0 + h*k3, ord))

  z1 = z0 + h*(k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine

subroutine runge_kutta_new_vector(tvec1,z0,tvec0,h)
  implicit none
  real(DP), intent(in) :: h
  complex(DP), intent(in) :: z0, tvec0
  complex(DP), intent(out) :: tvec1
  complex(DP) :: k1, k2, k3, k4, z1
  integer,parameter :: ord = 2

  ! k1 = conjg(tvec0 * action_val(z0, ord))
  ! k2 = conjg(tvec0 * action_val(z0 + 0.5_DP*h*k1, ord))
  ! k3 = conjg(tvec0 * action_val(z0 + 0.5_DP*h*k2, ord))
  ! k4 = conjg(tvec0 * action_val(z0 + h*k3, ord))

  k1 = conjg(action_val(z0, ord) * tvec0)
  call runge_kutta_new_config(z1, z0, 0.5_DP*h)
  k2 = conjg(action_val(z1, ord) * (tvec0 + 0.5_DP*h*k1))
  k3 = conjg(action_val(z1, ord) * (tvec0 + 0.5_DP*h*k2))
  call runge_kutta_new_config(z1, z0, h)
  k4 = conjg(action_val(z1, ord) * (tvec0 + h*k3))

  tvec1 = tvec0 + h*(k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine
end module

program generate_thimble
  use sub_generate_thimble
  implicit none
  type(critical_point) :: cp
  real(DP) :: t_init, dt
  real(DP) :: ori, w
  complex(DP) :: z0, z1
  complex(DP) :: tvec0, tvec1

  integer :: itre, ii

  ! relevant critical point
  cp%val = cmplx(-0.97492178167965704_DP, -0.53953762329131794_DP)
  cp%vec = cmplx(0.9603569245916534_DP, -0.2787733441146431_DP) ! thimble version
  ! cp%vec = cmplx(0.2787733514816438_DP, 0.9603569224531533_DP) ! anti-thimble version
  cp%kap = 1.9647512804566389_DP

  ! unrelevant critical point 1
  ! cp%val = cmplx(0.2777396340910929_DP, 2.1508620391644904_DP)
  ! cp%vec = cmplx(0.1616495103675127_DP, -0.9868482334168428_DP) ! thimble version
  ! cp%vec = cmplx(0.9868482334168428_DP, 0.1616495103675127_DP) ! anti-thimble version
  ! cp%kap = 3.7447743397899673_DP

  ! unrelevant critical point 2
  ! cp%val = cmplx(0.6971821475885641_DP, -1.6113244158731723_DP)
  ! cp%vec = cmplx(0.5277156596406994_DP, 0.8494210867231761_DP) ! thimble version
  ! cp%vec = cmplx(0.8494210867231761_DP, -0.5277156596406994_DP) ! anti-thimble version
  ! cp%kap = 2.5061451795334415_DP

  t_init = -10.0_DP
  ori = 1.0_DP
  call set_init_cond(z0, tvec0, cp, t_init, ori)

  itre = 1000000
  dt = 0.001_DP

  z1 = cmplx(0.0_DP, 0.0_DP)
  tvec1 = cmplx(0.0_DP, 0.0_DP)

  write(*,"('Runge Kutta time: ',ES24.15)") dt
  write(*,"('Init time: ', ES24.15, ' Init orientation', ES24.15)") t_init, ori
  write(*,"('Init config: ', 2ES24.15)") z0

  do ii = 1, itre
    !
    ! runge kuatta method (4th order)
    !
    call runge_kutta_new_config(z1, z0, dt)
    call runge_kutta_new_vector(tvec1, z0, tvec0, dt)

    write(*,'("time_config_tvec", 20ES24.15)') t_init+dt*ii, z1, tvec1, tvec1/abs(tvec1)
    write(*,'("error_check",10ES24.15)') &
      &aimag(action_val(cp%val, 0)) - aimag(action_val(z1, 0)), conjg(action_val(z1, 1)) - tvec1*cp%kap*ori

    !
    ! update variable
    !
    z0 = z1
    tvec0 = tvec1

    if (abs(z0) > 10.0_DP) stop
    if (abs(tvec0) > 10.0_DP) stop
  enddo

  stop
end program
