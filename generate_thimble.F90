module sub_generate_thimble
  implicit none
  integer, parameter :: DP=kind(1.0d0)
  real(DP), parameter :: PI=3.14159265359_DP
  complex(DP), parameter :: zi=(0.0_DP,1.0_DP)
  complex(DP), parameter :: sig=(1.0_DP, 0.0_DP)
  complex(DP), parameter :: lam=(0.333333333333333_DP, 0.0_DP)
  complex(DP), parameter :: exh=(1.0_DP, 1.0_DP)

contains

subroutine make_K_matrix(K,z)
  implicit none
  complex(DP),intent(inout) :: K(:,:)
  complex(DP),intent(in) :: z

  K(:,:) = sig + 3.0_DP*lam*z**2

  return
end subroutine

subroutine calc_diag_matrix(K, kap)
  use forpy_mod
  implicit none
  complex(DP),intent(in) :: K(:,:)
  real(DP),dimension(:,:),pointer,intent(out) :: kap

  integer :: ierr
  type(module_py) :: python_mod
  type(ndarray) :: nd_K, nd_kap
  type(tuple) :: arg
  type(list) :: paths
  type(object) :: ret

  ierr = forpy_initialize()
  ierr = get_sys_path(paths)
  ierr = paths%append(".")

  ierr = import_py(python_mod, "python_mod")

  ierr = ndarray_create(nd_K,K)
  ierr = tuple_create(arg,1)
  ierr = arg%setitem(0,nd_K)

  ierr = call_py(ret,python_mod,"takagi_factorization_diag",arg)
  ierr = cast(nd_kap,ret)
  ierr = nd_kap%get_data(kap)

  call python_mod%destroy
  call nd_K%destroy
  call nd_kap%destroy
  call arg%destroy
  call paths%destroy
  call ret%destroy
  call forpy_finalize()

  return
end subroutine

subroutine calc_ortho_vec(K,vec)
  use forpy_mod
  implicit none
  complex(DP),intent(in) :: K(:,:)
  complex(DP),dimension(:,:),pointer,intent(out) :: vec

  integer :: ierr
  type(module_py) :: python_mod
  type(ndarray) :: nd_K, nd_vec
  type(tuple) :: arg
  type(list) :: paths
  type(object) :: ret

  ierr = forpy_initialize()
  ierr = get_sys_path(paths)
  ierr = paths%append(".")

  ierr = import_py(python_mod, "python_mod")

  ierr = ndarray_create(nd_K,K)
  ierr = tuple_create(arg,1)
  ierr = arg%setitem(0,nd_K)

  ierr = call_py(ret,python_mod,"takagi_factorization_ortho_vec",arg)
  ierr = cast(nd_vec,ret)
  ierr = nd_vec%get_data(vec)

  call python_mod%destroy
  call nd_K%destroy
  call nd_vec%destroy
  call arg%destroy
  call paths%destroy
  call ret%destroy
  call forpy_finalize()

  return
end subroutine

subroutine set_init_cond(z, tvec, z_sig, vec, kap, uni, t_init)
  implicit none
  complex(DP), intent(in) :: z_sig,vec(:)
  real(DP), intent(in) :: uni(:),kap(:,:),t_init
  complex(DP), intent(inout) :: z,tvec(:)
  complex(DP) :: sum_grad
  integer :: i

  sum_grad = cmplx(0.0_DP, 0.0_DP)
  do i = 1, 1
    sum_grad = sum_grad + vec(i) * exp(kap(i,i)*t_init) * uni(i)
    tvec(i) = vec(i) * exp(kap(i,i)*t_init)
  end do

  z = z_sig + sum_grad

  return
end subroutine

function action_val(z,ord) result (S)
  implicit none
  complex(DP),intent(in) :: z
  integer,intent(in) :: ord
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

subroutine runge_kutta_new_config(z1,z0,h)
  implicit none
  real(DP),intent(in) :: h
  complex(DP),intent(in) :: z0
  complex(DP),intent(out) :: z1
  complex(DP) :: k1,k2,k3,k4
  integer,parameter :: ord=1

  k1 = h * conjg(action_val(z0, ord))
  k2 = h * conjg(action_val(z0 + 0.5_DP*k1, ord))
  k3 = h * conjg(action_val(z0 + 0.5_DP*k2, ord))
  k4 = h * conjg(action_val(z0 + k3, ord))

  z1 = z0 + (k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine

subroutine runge_kutta_new_vector(tvec1,z0,tvec0,h)
  implicit none
  real(DP),intent(in) :: h
  complex(DP),intent(in) :: z0,tvec0(:)
  complex(DP),intent(out) :: tvec1(:)
  complex(DP) :: k1,k2,k3,k4
  integer,parameter :: ord=2

  k1 = h * conjg(tvec0(1) * action_val(z0,ord))
  k2 = h * conjg(tvec0(1) * action_val(z0 + 0.5_DP*k1, ord))
  k3 = h * conjg(tvec0(1) * action_val(z0 + 0.5_DP*k2, ord))
  k4 = h * conjg(tvec0(1) * action_val(z0 + k3, ord))

  tvec1(:) = tvec0(:) + (k1 + 2.0_DP*k2 + 2.0_DP*k3 + k4)/6.0_DP

  return
end subroutine

function weight(t,z) result(w)
  implicit none
  real(DP),intent(in) :: t
  complex(DP),intent(in) :: z
  real(DP) :: w

  w = exp(-action_val(z,1))

  return
end function
end module

program generate_thimble
  use sub_generate_thimble
  implicit none
  real(DP) :: t
  complex(DP) :: z0,z1
  real(DP),dimension(:,:),pointer :: kap
  complex(DP) :: vec(1),tvec0(1),tvec1(1)
  real(DP) :: uni(1), w
  complex(DP) :: z_sig,K_mat(1,1)

  integer :: itre
  real(DP) :: h

  integer :: ii,iout

  z_sig = cmplx(-0.9749217865146699_DP, -0.5395376217014486_DP)

  !
  ! make matrix (K_mat) and calculate diagonal matrix (kappa) and orthonormal vector (vec)
  !
  call make_K_matrix(K_mat,z_sig)
  call calc_diag_matrix(K_mat,kap)
  ! call calc_ortho_vec(K_mat,vec)
  ! write(*,*) vec

  uni(1) = 1.0_DP
  vec(1) = cmplx(9.6035692459165345E-01, -2.7877334411464305E-01)

  t = -5.0_DP
  call set_init_cond(z0, tvec0, z_sig, vec, kap, uni, t)

  itre = 1000000
  h = 0.00001_DP
  ! itre = 100
  ! h = 0.1_DP

  iout = 99
  open(iout,file='thimble_generated_by_runge_kutta.dat',status='replace',form='formatted')
  write(iout,'(3ES24.15)') t,z0,tvec0(1)

  do ii = 1, itre
    !
    ! runge kuatta method (4th order)
    !
    call runge_kutta_new_config(z1,z0,h)
    call runge_kutta_new_vector(tvec1,z0,tvec0,h)

    write(*,'(5ES24.15)') t,z1,tvec1(1)
    ! write(*,'(3ES24.15)') t, conjg(action_val(z1,1)) - tvec1(1)*kap(1,1)

    !
    ! update variable
    !
    t = t + h
    z0 = z1
    tvec0(:) = tvec1(:)

    write(iout,'(5ES24.15)') t,z0,tvec0(1)

    w = weight(t,z0)
    if (w <= 1E-10) stop 'weight is small than 10^(-10)'

  enddo
  close(iout)

  stop
end program
