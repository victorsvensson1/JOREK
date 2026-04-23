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

    subroutine initialise_on_flux_surfaces(sim, R_mat, Z_mat)
        type(particle_sim), intent(inout) :: sim
        real(kind=8), intent(in)          :: R_mat(:,:), Z_mat(:,:)
        integer :: s, idx, n_psi, ifail
        real(kind=8) :: r_start, z_start
        real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt, B_star
        real(kind=8) :: normB

        ! We only care about the number of radial surfaces (n_psi)
        n_psi = size(R_mat, 2)

        select type (p => sim%groups(1)%particles)
        type is (particle_gc_relativistic)
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
                p(idx)%q = -1
                p(idx)%i_elm = -1 


                call find_RZ(sim%fields%node_list, sim%fields%element_list, &
                                p(idx)%x(1), p(idx)%x(2), &      ! Input R, Z for particle idx
                                p(idx)%x(1), p(idx)%x(2), &      ! Dummy/Initial R, Z
                                p(idx)%i_elm, &                  ! The specific particle's element
                                p(idx)%st(1), p(idx)%st(2), &    ! The specific particle's local coords
                                ifail)
                
                call sim%fields%calc_EBNormBGradBCurlbDbdt(sim%time,p(idx)%i_elm,p(idx)%st, p(idx)%x(3),E,b,normB,gradB,curlb,dbdt)
                
                print *, "E: ", E
                print *, "dot_prod:", dot_product(E, b)

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
        
        real(kind=8)    :: local_time, t_step
        integer(kind=4) :: j, k, n_pushes, n_lost

        n_lost = 0 
        t_step = timesteps(1)
        
        ! Calculate total pushes required to reach end_time
        ! Or you can set n_pushes to a fixed integer (e.g., 1000)
        n_pushes = ceiling((end_time - sim%time) / t_step)
        print *, "Running particle trace for ", n_pushes, " pushes per particle to reach end time: ", end_time

        select type (particles => sim%groups(1)%particles)
        type is (particle_gc_relativistic)

            !!$omp parallel do default(private) shared(sim, t_step, n_pushes) reduction(+:n_lost)
            do j = 1, size(particles, 1)
                if (particles(j)%i_elm <= 0) cycle 

                local_time = sim%time

                ! Inner Loop: Perform exactly N pushes for particle j
                do k = 1, n_pushes
                    
                    call runge_kutta_fixed_dt_gc_push_jorek(sim%fields, local_time, t_step, sim%groups(1)%mass, particles(j))
                    
                    local_time = local_time + t_step

                    ! If particle is lost, stop pushing it and move to next particle (j+1)
                    if (particles(j)%i_elm <= 0) then
                        n_lost = n_lost + 1
                        exit 
                    end if
                end do
            end do
            !!$omp end parallel do

        end select

        ! Update global simulation time based on the steps taken
        !sim%time = sim%time + (n_pushes * t_step)

        print *, "Tracing complete."
        print *, "Pushes per particle:", n_pushes
        print *, "Total particles lost:", n_lost
        print *, "Current simulation time:", sim%time
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
