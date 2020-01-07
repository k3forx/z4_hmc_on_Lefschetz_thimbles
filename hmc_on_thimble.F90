module hmc_on_thimble_mod
  use generate_thimble_mod
  implicit none

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
  real(DP), intent(out) :: t_init, ori
  real(DP) :: grand

  t_init = -5.0_DP

  grand = generate_gauss_rand()
  ori = grand / sqrt(grand**2)

  return
end subroutine

subroutine calc_inverse_matrix(mat_sca, mat_inv_sca)
  implicit none
  complex(DP), intent(in) :: mat_sca
  complex(DP), intent(inout) :: mat_inv_sca
  complex(DP) :: mat(1,1), mat_inv(1,1)
  integer, allocatable :: IPIV(:)
  complex(DP), allocatable :: WORK(:)
  integer :: N, LDA, LWORK, INFO

  mat(1,1) = mat_sca
  mat_inv(1,1) = mat(1,1)

  N = int(sqrt(real(size(mat), kind=DP)))
  LDA = max(1, N)
  LWORK = max(1, N)

  allocate (IPIV(N), WORK(max(1, LWORK)))

  call ZGETRF(N, N, mat_inv, N, IPIV, INFO)
  if (INFO /= 0) stop 'ZGETRF is failed'

  call ZGETRI(N, mat_inv, LDA, IPIV, WORK, LWORK, INFO)
  if (INFO /= 0) stop 'ZGETRI is failed'

  mat_inv_sca = mat_inv(1,1)
  deallocate (IPIV, WORK)
  return
end subroutine

subroutine generate_gaussian_momentum(p, tvec, inv)
  implicit none
  complex(DP), intent(inout) :: p
  complex(DP), intent(in) :: tvec, inv

  p = cmplx(generate_gauss_rand(), generate_gauss_rand())
  p = tvec*real(inv*p, kind=DP)

  return
end subroutine

function calcu_hamil(z, p) result (hamil)
  implicit none
  complex(DP), intent(in) :: z, p
  real(DP) :: hamil

  hamil = 0.5_DP*(action_val(z, 0) + conjg(action_val(z, 0)))
  hamil = hamil + 0.5_DP*p*conjg(p)

  return
end function

subroutine find_next_candidate(dtau, t_init, ori0, t0, z0, p0, tvec0, inv0, cp, ori1, t1, z1, p1, tvec1, inv1)
  implicit none
  type(critical_point), intent(in) :: cp
  real(DP), intent(in) :: dtau, t_init, t0, ori0
  complex(DP), intent(in) :: z0, p0, tvec0, inv0
  real(DP), intent(inout) :: t1, ori1
  complex(DP), intent(inout) :: z1, p1, tvec1, inv1

  real(DP) :: prod, ori0_tmp, t0_tmp
  complex(DP) :: z_tmp, tvec_tmp
  real(DP) :: error, sum_(2), consts(2)
  integer :: ndim, ii, jj, iter

  error = 1.0_DP
  t0_tmp = t0
  ori1 = ori0

  do while (error > 1E-10)
    !
    ! set initial condition z_i(t_0) and solve flow eq. to calculate z_i[e_k, t'_k]
    !
    call set_init_cond(z_tmp, tvec_tmp, cp, t_init, ori0)
    call solve_flow_eq(z_tmp, z1, tvec_tmp, tvec1, t0_tmp)
    call calc_inverse_matrix(tvec1, inv1)
    ! write(*,'(10ES24.15)') t0_tmp, z0, z1, z0-z1

    write(*,'("before ",20ES24.15)') ori0, t0_tmp, error, z1, p1, tvec1, inv1

    !
    ! calculate candidate of time t'_(k+1)
    !
    prod = real(inv0*(z0 + dtau*p0 - 0.5_DP*dtau**2*conjg(action_val(z0, 1)) - z1), kind=DP)
    t1 = t0_tmp + prod / cp%kap*ori0
    if (t1 < 0.0_DP) exit

    !
    ! evaluate error
    !
    p1 = tvec0*ori0*cp%kap*(t1 - t0_tmp)
    error = p1*conjg(p1)

    !
    ! update orientation vector and time
    !
    t0_tmp = t1

    write(*,'("after  ",20ES24.15)') ori0, t0_tmp, error, z1, p1, tvec1, inv1
  end do

  !
  ! grid search method
  !
  if (t1 < 0.0_DP) then
    write(*,'("ANOTHER ITERATION STARTS !!!")')
    call itrerate_calc_candidate_time(t1, z0, z1, p0, tvec0, tvec1, inv0, dtau, cp, consts(1), ori0)
    call set_init_cond(z_tmp, tvec_tmp, cp, t_init, ori0)
    call solve_flow_eq(z_tmp, z1, tvec_tmp, tvec1, t1)
    call calc_inverse_matrix(tvec1, inv1)

    consts(1) = 2.0_DP*aimag(inv0*(z0 - z1))/dtau**2
    consts(2) = 2.0_DP*aimag(inv1*(p1 - 0.5_DP*dtau*conjg(action_val(z1, 1))))/dtau
    p1 = p1 - 0.5_DP*dtau*conjg(action_val(z1, 1)) - 0.5_DP*dtau*zi*(tvec1*consts(2))

    ! write(*,'(20ES24.15)') z1, z0 + dtau*p0 - 0.5_DP*dtau**2*conjg(action_val(z0, 1)) - 0.5_DP*dtau**2*zi*tvec0*consts(1)
    ! write(*,'("after  ",10ES24.15)') ori0_tmp, t0_tmp, error, z1, tvec1, p1
    write(*,'("config_and_force", 40ES24.15)') &
      &z0, dtau*p0 - 0.5_DP*dtau**2*conjg(action_val(z0, 1)), -0.5_DP*dtau**2*zi*tvec0*consts(1), z1, t1
    write(*,'("error ",10ES24.15)') &
      &z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)) + 0.5_DP*dtau**2*zi*tvec0*consts(1)
    return
  end if

  call calc_inverse_matrix(tvec1, inv1)

  !
  ! calculate one of Lagrange multiplies lamda[r], w^(n+1/2) and z^(n+1)
  !
  consts(1) = 2.0_DP*aimag(inv0*(z0 - z1))/dtau**2
  p1 = p0 - 0.5_DP*dtau*conjg(action_val(z0, 1)) - 0.5_DP*dtau*zi*(tvec0*consts(1))
  z1 = z0 + dtau*p1

  !
  ! calculate inverse matrix of tangent vector and another Langrange multiplier lambda[v] and w^(n+1)
  !
  consts(2) = 2.0_DP*aimag(inv1*(p1 - 0.5_DP*dtau*conjg(action_val(z1, 1))))/dtau
  p1 = p1 - 0.5_DP*dtau*conjg(action_val(z1, 1)) - 0.5_DP*dtau*zi*(tvec1*consts(2))

  write(*,'("config_and_force", 40ES24.15)') &
    &z0, dtau*p0 - 0.5_DP*dtau**2*conjg(action_val(z0, 1)), -0.5_DP*dtau**2*zi*tvec0*consts(1),&
    &(dtau*p0 - 0.5_DP*dtau**2*conjg(action_val(z0, 1)))/tvec0, tvec0
  write(*,'("error ",10ES24.15)') &
    &z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)) + 0.5_DP*dtau**2*zi*tvec0*consts(1)

  return
end subroutine

subroutine itrerate_calc_candidate_time(t_cand, z0, z1, p0, tvec0, tvec1, inv0, dtau, cp, lamr, ori)
  implicit none
  complex(DP), intent(in) :: z0, p0, tvec0, inv0
  complex(DP), intent(inout) :: z1, tvec1
  real(DP), intent(in) :: dtau, ori
  real(DP), intent(inout) :: lamr, t_cand
  type(critical_point), intent(in) :: cp
  complex(DP) :: z_init, tvec_init
  real(DP) :: error, dt, t1, d0, d1, t_store
  real(DP), parameter :: t_init = -5.0_DP
  integer :: iter

  error = abs(z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)))
  dt = 0.01_DP
  t_cand = t_init
  lamr = 0.0_DP

  d0 = abs(- dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)))

  do while (abs(d1 - d0) >= 1E-14)
    call calc_candidate_time(t_cand, t1, dt, lamr, dtau, z0, p0, tvec0, inv0, cp, ori)
    call set_init_cond(z_init, tvec_init, cp, t_init, ori)
    call solve_flow_eq(z_init, z1, tvec_init, tvec1, t_cand)

    lamr = 2.0_DP*aimag(inv0*(z0 - z1))/dtau**2
    d1 = abs(z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)) + 0.5_DP*dtau**2*zi*tvec0*lamr)
    write(*,'("iter end", I10, 10ES24.15)') iter, t_cand, lamr, z1, d0, d1

    if (abs(d1 - d0) < 1E-14) exit

    t_store = t_cand
    t_cand = t_cand - dt
    dt = dt * 0.1_DP
    d0 = d1

  end do

  t_cand = t_store

  return
end subroutine

subroutine calc_candidate_time(t0, t1, dt, lamr, dtau, z0, p0, tvec0, inv0, cp, ori)
  implicit none
  real(DP), intent(inout) :: t0, t1, dt
  real(DP), intent(in) :: lamr, dtau, ori
  complex(DP), intent(in) :: z0, p0, tvec0, inv0
  type(critical_point), intent(in) :: cp

  real(DP) :: t_init, t_store
  real(DP) :: d0, d1
  integer :: iter
  complex(DP) :: z_init, tvec_init, z1, tvec1

  t_init = -5.0_DP

  d0 = 1.0_DP
  d1 = 0.0_DP
  z1 = z0

  iter = 0
  do while (d1 - d0 <= 0.0_DP)
    call set_init_cond(z_init, tvec_init, cp, t_init, ori)
    call solve_flow_eq(z_init, z1, tvec_init, tvec1, t0)
    d0 = abs(z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)) + 0.5_DP*dtau**2*zi*tvec0*lamr)

    t1 = t0 + dt
    call set_init_cond(z_init, tvec_init, cp, t_init, ori)
    call solve_flow_eq(z_init, z1, tvec_init, tvec1, t1)
    d1 = abs(z1 - z0 - dtau*p0 + 0.5_DP*dtau**2*conjg(action_val(z0, 1)) + 0.5_DP*dtau**2*zi*tvec0*lamr)

    write(*,'("iter ",I10,10ES24.15)') iter, t0, t1, z1, lamr, d1 - d0
    t_store = t0
    t0 = t1
    iter = iter + 1

    if ((iter > 1E+6) .or. ((d1 - d0) == 0.0_DP)) exit
  end do

  t0 = t_store

  return
end subroutine
end module

program hmc_on_thimble
  use mt19937
  use hmc_on_thimble_mod
  implicit none
  type(critical_point) :: cp
  real(DP) :: ori0, ori1, ori2, ori3, ori0_tmp
  real(DP) :: t_init, t0, t1, t2, t3, t0_tmp
  complex(DP) :: z0, z1, z2, z3, z0_tmp
  complex(DP) :: tvec0, tvec1, tvec2, tvec3, tvec0_tmp
  complex(DP) :: inv0, inv1, inv2, inv3, inv0_tmp
  complex(DP) :: p0, p1, p2, p3, p0_tmp
  real(DP) :: rho, h0, h1, h2, rand, pacc
  integer :: ndim, ii

  integer :: NTHERM, NSKIP, NSAMPLE, nmd, seed
  real(DP) :: dtau

  integer :: istep, jstep, iout, itry, iacc
  real(DP) :: dist_from_cp, sum_of_p_and_f, costheta
  real(DP) :: lamr

  iout = 99
  open(iout, file='HMC_PARAM', status='old', form='formatted', action='read')
  read(iout,*) NTHERM, NSKIP, NSAMPLE
  read(iout,*) nmd, seed
  close(iout)

  call sgrnd(seed)

  cp%val = cmplx(-0.97492178167965704_DP, -0.53953762329131794_DP)
  cp%vec = cmplx(0.9603569245916534_DP, -0.2787733441146431_DP)
  cp%kap = 1.9647512804566389_DP

  dtau = 1.0_DP / real(nmd, kind=DP)

  !
  ! initialization procedure
  !
  call set_init_param(t_init, ori0)
  call set_init_cond(z0, tvec0, cp, t_init, ori0)

  !
  ! solve flow equation via Runge Kutta method
  !
  t0 = -t_init
  call solve_flow_eq(z0, z1, tvec0, tvec1, t0)

  itry = 0
  iacc = 0
  do istep = 1, NTHERM + NSKIP*NSAMPLE
    !
    ! calculate inverse matrix and generate initial momentum
    !
    call calc_inverse_matrix(tvec0, inv0)
    call generate_gaussian_momentum(p0, tvec0, inv0)

    ! nmd = 10, before fliped
    ! ori0 = 1.000000000000000_DP
    ! t0 = 4.100374401268285_DP
    ! z0 = cmplx(-0.8248587560409677_DP, -0.5846964407533598_DP)
    ! p0 = cmplx(-1.640588982374780_DP, 0.5109663777076791_DP)
    ! tvec0 = cmplx(0.1374207201259041_DP, -0.04280003336325399_DP)
    ! inv0 = cmplx(6.633322675583814_DP, 2.065968317550841_DP)

    ! nmd = 10, after fliped
    ori0 = -1.000000000000000E+00
    t0 = 2.908322304880164_DP
    z0 = cmplx(-0.9908249224939001_DP, -0.5349393609008363_DP)
    p0 = cmplx(-1.663213881030954_DP, 0.4790086761836458_DP)
    tvec0 = cmplx(-0.01605274907526825_DP, 0.004623223878681929_DP)
    inv0 = cmplx(-57.51949068446022_DP, -16.56571978009368_DP)
    p0 = -p0

    ! nmd = 10, after filped 2
    ! ori0 = -1.000000000000000_DP
    ! t0 = 4.446329999999882_DP
    ! z0 = cmplx(-1.375754293881792_DP, -0.4346375542349374_DP)
    ! p0 = cmplx(-0.3539262112351909_DP, 0.08270528152578302_DP)
    ! tvec0 = cmplx(-0.4987677566642592_DP, 0.1165517738483881_DP)
    ! inv0 = cmplx(-1.901127995314834_DP,  -0.4442545397254464_DP)

    !
    ! Compute initial hamiltonian
    !
    h0 = calcu_hamil(z0, p0)

    !
    ! calculate z^(n+1) and w^(n+1)
    !
    do jstep = 1, nmd
      call set_init_cond(z0, tvec0, cp, t_init, ori0)
      call solve_flow_eq(z0, z1, tvec0, tvec1, t0)
      write(*,'("config_and_moment",20ES24.15)') ori0, t0, z0, p0, tvec0, inv0

      !
      ! calculate distance from critical point
      !
      p1 = p0 - 0.5_DP*dtau*conjg(action_val(z0, 1))
      z1 = z0 + dtau*p1
      sum_of_p_and_f = abs(dtau*p1)
      costheta = real(abs(dtau*p1)*(cp%val - z0)/(abs(cp%val - z0)*dtau*p1), kind=DP)

      dist_from_cp = abs(z0 - cp%val)

      if (0 < costheta .and. costheta < 1) then
        if (dist_from_cp < sum_of_p_and_f) then
          t0 = t0 - log(dist_from_cp / (abs(z1 - z0) - dist_from_cp)) / cp%kap
          ori0 = - ori0
          write(*,'("FLIPED !!!!")')
          write(*,'("Fliped_orientatin: ",ES24.15," Approximation_time: ",ES24.15)') ori0, t0
          ! stop
        end if
      end if

      if (jstep == 2) stop
      call find_next_candidate(dtau, t_init, ori0, t0, z0, p0, tvec0, inv0, cp, ori1, t1, z1, p1, tvec1, inv1)

      ori0 = ori1
      t0 = t1
      z0 = z1
      p0 = p1
      tvec0 = tvec1
      inv0 = inv1
    end do
    write(*,'("config_and_moment",20ES24.15)') ori0, t0, z0, p0, tvec0, inv0

#ifdef _CHECK_REVERSE
    ori2 = ori1
    t2 = t1
    z2 = z1
    p2 = -p1
    tvec2 = tvec1
    inv2 = inv1

    do jstep = 1, nmd
      write(*,'("@config_and_moment",20ES24.15)') ori2, t2, z2, p2, tvec2, inv2
      !
      ! calculate distance from critical point
      !
      ! p1 = p0 - 0.5_DP*dtau*conjg(action_val(z0, 1))
      ! z1 = z0 + dtau*p1
      ! sum_of_p_and_f = abs(dtau*p1)
      ! costheta = real(abs(dtau*p1)*(cp%val - z0)/(abs(cp%val - z0)*dtau*p1), kind=DP)

      ! dist_from_cp = abs(z0 - cp%val)

      ! if (0 < costheta .and. costheta < 1) then
      !   if (dist_from_cp < sum_of_p_and_f) then
      !     t0 = t0 - log(dist_from_cp / (abs(z1 - z0) - dist_from_cp)) / cp%kap
      !     ori0 = - ori0
      !     write(*,'("FLIPED !!!!")')
      !     stop
      !   end if
      ! end if

      call find_next_candidate(dtau, t_init, ori2, t2, z2, p2, tvec2, inv2, cp, ori3, t3, z3, p3, tvec3, inv3)

      ori2 = ori3
      t2 = t3
      z2 = z3
      p2 = p3
      tvec2 = tvec3
      inv2 = inv3
    end do

    h2 = calcu_hamil(z1, p1)
    write(*,'("@config_and_moment",20ES24.15)') ori2, t2, z2, p2, tvec2, inv2
    write(*,'(20ES24.15)') z0 - z2, p0 - p2
    write(*,'("hamil", 10ES24.15)') (1.0_DP/real(nmd, kind=DP))**2, h0, h2, h0 - h2

    stop
#endif

    h1 = calcu_hamil(z1, p1)

    rho = min(1.0_DP, exp(h0-h1))
    rand = grnd()

    itry = itry + 1
    if (rand <= rho) then

      iacc = iacc + 1
      continue

    else

      z1 = z0
      ori1 = ori0
      t1 = t0
      tvec0 = tvec1

    end if

    ! write(*,'("hmc_config_and_moment",21ES24.15)') ori1, t1, z1, p1, tvec1, inv1

    ori0 = ori1
    t0 = t1
    z0 = z1
    tvec0 = tvec1

    if (ori0 < 0.0_DP) stop
  end do


  !
  ! Out put result of HMC simulation
  !
  pacc = real(iacc, kind=DP)/itry
  write(*,'("# HMC Metropolis test statistics.")')
  write(*,'("# dt= ",ES14.6," itry =",I10," iacc =",I10," Pacc = ",F10.6)') &
      & dtau,itry,iacc,pacc*100

  stop
end program



