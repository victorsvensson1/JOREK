!> Module for particle tracing in JOREK to keep track of flux surfaces.
module mod_fluxsurf_evol
    use particle_tracer
    ! Use the standard (non-relativistic) guiding center module
    use mod_find_rz_nearby
    use mod_gc_relativistic
    use equil_info
    use mod_fields_linear  
    use mod_expression    
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use settings

    implicit none

    private
    public :: initialise_on_flux_surfaces, run_particle_trace, get_psi_at_pos

contains

    subroutine initialise_on_flux_surfaces(sim, R_mat, Z_mat)
        type(particle_sim), intent(inout) :: sim
        real(kind=8), intent(in)          :: R_mat(:,:), Z_mat(:,:)
        integer :: s, idx, n_psi, ifail
        real(kind=8) :: r_start, z_start

        ! We only care about the number of radial surfaces (n_psi)
        n_psi = size(R_mat, 2)

        select type (p => sim%groups(1)%particles)
        type is (particle_fieldline)
            idx = 0
            ! Loop only over flux surfaces
            do s = 1, n_psi
                idx = idx + 1
                
                ! Use the first theta entry (theta = 0) for each surface
                r_start = R_mat(1, s)
                z_start = Z_mat(1, s)
                
                p(idx)%x(1) = r_start
                p(idx)%x(2) = z_start
                p(idx)%x(3) = 0.0d0

                p(idx)%v = 1.0d5

                ! Initialize element search
                p(idx)%i_elm = -1 
                call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, &
                    r_start, z_start, 0.5d0, 0.5d0, p(idx)%i_elm, &
                    r_start, z_start, p(idx)%st(1), p(idx)%st(2), p(idx)%i_elm, ifail)
                
                if (ifail /= 0) then
                    p(idx)%i_elm = -1
                    print *, "Warning: Particle at surface ", s, " failed localization."
                endif
            end do
            print *, "Initialized ", idx, " particles (one per flux surface at theta=0)."
        end select
    end subroutine initialise_on_flux_surfaces

    subroutine run_particle_trace(sim, timesteps, events, end_time)
        type(particle_sim), intent(inout) :: sim
        real(kind=8), intent(in)          :: timesteps(:), end_time
        type(event), intent(inout)        :: events(:)
        
        real(kind=8)    :: target_time, t_step, current_time
        integer(kind=4) :: i, j, k, n_steps, n_lost

        ! Process any events at the very start (t=0)
        !call with(sim, events, at=sim%time)
        current_time = sim%time

        do while (current_time < end_time)
            
            target_time = current_time + minval(timesteps) 
            if (target_time > end_time) target_time = end_time

            !do i = 1, size(sim%groups)
            !n_steps = nint((target_time - current_time) / timesteps(i))
            n_steps = 1
            if (n_steps <= 0) cycle
            
            n_lost = 0  
            select type (particles => sim%groups(1)%particles)
            type is (particle_fieldline)
        
                !$omp parallel do default(private) shared(sim, n_steps, timesteps, i, target_time) reduction(+:n_lost)
                do j = 1, size(particles, 1)
                    if (particles(j)%i_elm <= 0) cycle 
                    do k = 1, n_steps
                        t_step = current_time + (k-1)*timesteps(i)
                        !call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, particles(j), t_step, timesteps(i))
                        call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, particles(j), t_step, timesteps(1))
                        if (particles(j)%i_elm <= 0) exit
                    end do 
                end do
                !$omp end parallel do
            end select
            !end do

            ! Move the simulation clock forward
            current_time = target_time
            
        end do

        print *, "Tracing complete. Returning to main script at time: ", sim%time
    end subroutine run_particle_trace

    subroutine get_psi_at_pos(ES, node_list, element_list, x, psi, ierr)

    ! --- Routine parameters
    type(t_equil_state), intent(in)      :: ES
    type(type_node_list), intent(in)     :: node_list
    type(type_element_list), intent(in)  :: element_list
    real(kind=8), intent(in)             :: x(3)   ! R, Z, Phi
    real(kind=8), intent(out)            :: psi
    integer, intent(out)                 :: ierr

    ! --- Local variables
    integer :: units
    type(t_expr_list) :: expr_list
    real*8, allocatable, save :: result(:,:,:,:), res0d(:)
    character(len=16)  :: name_psi
    character(len=128) :: desc_psi

    ierr = 0

    name_psi = 'Psi'
    desc_psi = 'poloidal magnetic flux'

    call add(expr_list, name_psi, desc_psi)

    units = get_int_setting('units', ierr)

    call eval_expr(ES, units, expr_list, &
                   pol_pos(node_list, element_list, ES, R=x(1), Z=x(2)), &
                   tor_pos(phi=x(3)), result, ierr)

    if (ierr /= 0) return

    call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)

    if (allocated(res0d)) then
        psi = res0d(1)
    else
        ierr = -1 ! Evaluation failed to produce result
    end if

end subroutine get_psi_at_pos

end module mod_fluxsurf_evol
