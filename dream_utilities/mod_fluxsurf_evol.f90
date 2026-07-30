!> Module for particle tracing in JOREK to keep track of flux surfaces.
module mod_fluxsurf_evol
    use particle_tracer
    ! Use the standard (non-relativistic) guiding center module
    !use mod_find_rz_nearby
    !use find_RZ
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

    subroutine initialise_on_flux_surfaces(sim, particles, R_mat, Z_mat)
        type(particle_sim), intent(inout)              :: sim
        type(particle_gc_relativistic), intent(inout)  :: particles(:)
        real(kind=8), intent(in)                       :: R_mat(:,:), Z_mat(:,:)
        integer :: s, n_psi, ifail
        real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt
        real(kind=8) :: normB

        n_psi = size(R_mat, 2)

        do s = 1, n_psi
            particles(s)%x(1) = R_mat(1, s)
            particles(s)%x(2) = Z_mat(1, s)
            particles(s)%x(3) = 0.0d0
            particles(s)%q    = -1
            particles(s)%i_elm = -1

            call find_RZ(sim%fields%node_list, sim%fields%element_list, &
                        particles(s)%x(1), particles(s)%x(2), &
                        particles(s)%x(1), particles(s)%x(2), &
                        particles(s)%i_elm, &
                        particles(s)%st(1), particles(s)%st(2), ifail)

            call sim%fields%calc_EBNormBGradBCurlbDbdt(sim%time, particles(s)%i_elm, &
                        particles(s)%st, particles(s)%x(3), E, b, normB, gradB, curlb, dbdt)

            if (ifail /= 0) then
                particles(s)%i_elm = -1
                print *, "Warning: Particle at surface ", s, " failed localization."
            endif
        end do
        print *, "Initialized ", n_psi, " particles (one per flux surface at theta=0)."
    end subroutine initialise_on_flux_surfaces


    subroutine run_particle_trace(sim, particles, mass, timesteps, end_time)
        type(particle_sim), intent(inout)              :: sim
        type(particle_gc_relativistic), intent(inout)  :: particles(:)
        real(kind=8), intent(in)                       :: mass
        real(kind=8), intent(in)                       :: timesteps(:), end_time

        real(kind=8)    :: local_time, t_step
        integer(kind=4) :: j, k, n_pushes, n_lost

        n_lost  = 0
        t_step  = timesteps(1)
        n_pushes = ceiling((end_time - sim%time) / t_step)
        print *, "Running particle trace for ", n_pushes, " pushes per particle"

        do j = 1, size(particles, 1)
            if (particles(j)%i_elm <= 0) cycle
            local_time = sim%time
            do k = 1, n_pushes
                call runge_kutta_fixed_dt_gc_push_jorek(sim%fields, local_time, t_step, mass, particles(j))
                local_time = local_time + t_step
                if (particles(j)%i_elm <= 0) then
                    n_lost = n_lost + 1
                    exit
                end if
            end do
        end do
        print *, "Tracing complete. Lost:", n_lost
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

    if (allocated(result)) deallocate(result)
    if (allocated(res0d))  deallocate(res0d)

    name_psi = 'Psi_N'
    desc_psi = 'normalized poloidal magnetic flux'

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
