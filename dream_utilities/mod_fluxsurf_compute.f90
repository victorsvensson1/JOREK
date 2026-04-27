!> Module to determine flux surfaces, inspired by the 'fluxsurface' command in exec_commands.f90.
!> Provides a subroutine to determine (R,Z) values of the flux surfaces to be given as input to DREAM

module mod_fluxsurf_compute
    use mod_fields_linear  ! For fields_linear type
    use mod_expression     ! For expression evaluation
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use equil_info
    use settings
    use convert_character
    implicit none

    !real*8, allocatable, private, save :: result(:,:,:,:), res0d(:)


    private
    public :: fluxsurface

contains

  subroutine fluxsurface(ES, node_list, element_list, psi_list, R_mat, Z_mat, theta_list, ierr)

    type(t_equil_state), intent(in)      :: ES
    type(type_node_list), intent(in)     :: node_list
    type(type_element_list), intent(in)  :: element_list
    real*8, intent(in)                   :: psi_list(:)
    integer, intent(out)                 :: ierr
    real*8, allocatable, intent(out)     :: R_mat(:,:)
    real*8, allocatable, intent(out)     :: Z_mat(:,:)
    real*8, allocatable, intent(out)     :: theta_list(:)

    integer :: i, k, iter, n_psi, n_theta
    real*8  :: target_psi, curr_theta, rho_0, rho_1, rho_next
    real*8  :: f0, f1, pi
    real*8, allocatable :: rho_prev(:)

    logical :: converged

    n_psi = size(psi_list)
    n_theta = 128 
    pi = 4.0d0 * atan(1.0d0)

    allocate(theta_list(n_theta), R_mat(n_theta, n_psi), Z_mat(n_theta, n_psi), rho_prev(n_theta))

    do i = 1, n_theta
        theta_list(i) = (real(i-1, 8) / real(n_theta-1, 8)) * 2.0d0 * pi
    end do

    do k = 1, n_psi
        target_psi = ES%psi_axis + psi_list(k) * (ES%psi_bnd - ES%psi_axis)

        do i = 1, n_theta
            curr_theta = theta_list(i)

            ! Set two initial guesses for the Secant Method
            if (k == 1) then
                rho_0 = 0.00d0
                rho_1 = 0.02d0
            else
                if (rho_prev(i) < 1.d-2) then
                    ! rho_prev collapsed to axis — use analytic estimate instead
                    ! psi ~ rho^2 near axis, so rho ~ sqrt(psiN) * a
                    rho_1 = sqrt(psi_list(k)) * ES%LCFS_a
                    rho_0 = rho_1 * 0.8d0

                else
                    rho_0 = rho_prev(i) * 0.98d0
                    rho_1 = rho_prev(i)
                end if
            end if


            ! Get initial psi values
            call get_psi_only(node_list, element_list, ES, &
                              ES%R_axis + rho_0*cos(curr_theta), &
                              ES%Z_axis + rho_0*sin(curr_theta), f0, ierr)

            f0 = f0 - target_psi

            converged = .false.

            do iter = 1, 20
                call get_psi_only(node_list, element_list, ES, &
                                  ES%R_axis + rho_1*cos(curr_theta), &
                                  ES%Z_axis + rho_1*sin(curr_theta), f1, ierr)
                f1 = f1 - target_psi

                if (abs(f1 - f0) < 1.d-15) exit ! Prevent division by zero

                ! Secant Method Formula
                rho_next = rho_1 - f1 * (rho_1 - rho_0) / (f1 - f0)
                rho_next = max(0.0d0, min(rho_next, sqrt(psi_list(k)) * ES%LCFS_a * 2.0d0))

                ! Update guesses
                rho_0 = rho_1
                f0 = f1
                rho_1 = rho_next

                if (abs(rho_1 - rho_0) < 1.d-8) then
                    converged = .true.
                    exit
                end if
            end do

            if (converged) then
                R_mat(i, k) = ES%R_axis + rho_1 * cos(curr_theta)
                Z_mat(i, k) = ES%Z_axis + rho_1 * sin(curr_theta)
                rho_prev(i) = rho_1
            else
                print *, "Warning: Secant method did not converge for k=", k, " theta index=", i, "theta=", curr_theta
                ! If it failed, use a linear interpolation between neighbors as a fallback
                rho_prev(i) = 0.5d0 * (rho_prev(max(1,i-1)) + rho_prev(min(n_theta,i+1)))
            end if
        end do
        print *, "DEBUG: Target Psi for k=", k, " is ", target_psi
        print *, "DEBUG: Resulting R at theta=0 is ", R_mat(1, k)
    end do

    deallocate(rho_prev)
  end subroutine fluxsurface

  subroutine get_psi_only(node_list, element_list, ES, R, Z, psi_val, ierr)
    type(type_node_list),    intent(in)  :: node_list
    type(type_element_list), intent(in)  :: element_list
    type(t_equil_state),     intent(in)  :: ES
    real*8,                  intent(in)  :: R, Z
    real*8,                  intent(out) :: psi_val
    integer,                 intent(out) :: ierr

    type(t_expr_list) :: expr_list_local
    integer :: units
    character(len=16)  :: name_psi
    character(len=128) :: desc_psi
    real*8, allocatable :: result(:,:,:,:), res0d(:)

    if (allocated(res0d)) deallocate(res0d)
    if (allocated(result)) deallocate(result)

    name_psi = 'Psi'
    desc_psi = 'Poloidal flux'

    ierr = 0
    units = get_int_setting('units', ierr)

    call add(expr_list_local, name_psi, desc_psi)

    call eval_expr(ES, units, expr_list_local, pol_pos(node_list, element_list, ES, R=R, Z=Z), tor_pos(nphi=4*n_plane), result, ierr)
    call apply_four_filter(result, simple_filter(n=0), expr_list_local%n_coord, ierr)
    call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)


    psi_val = res0d(1)

  end subroutine get_psi_only


end module mod_fluxsurf_compute
