module mod_com_dream
    use particle_tracer
    use mod_fluxsurf_evol  ! Assuming your tracer subroutines are here
    use mod_computeB       ! For comp_B_field
    use mod_fields_linear
    use equil_info
    use mod_fluxsurf_evol, only: initialise_on_flux_surfaces, run_particle_trace, get_psi_at_pos
    use mod_fluxsurf_avg, only: avg_fluxsurf_list
    use mod_fluxsurf_compute, only: fluxsurface
    use mod_computeB, only: comp_B_field
    use mod_save_flux_hdf5
    implicit none

    private 
    public :: com_dream

contains

    subroutine com_dream(sim, ES, n_psi_surfaces, dt_trace)
        type(particle_sim), intent(inout) :: sim
        type(t_equil_state), intent(in)   :: ES
        integer, intent(in)               :: n_psi_surfaces
        real(kind=8), intent(in)          :: dt_trace

        ! --- Local Data ---
        real(kind=8), allocatable :: PsiN_list(:), Psi_list(:), theta_list(:)
        real(kind=8), allocatable :: R_mat(:,:), Z_mat(:,:)
        real(kind=8), allocatable :: Br_mat(:,:), Bz_mat(:,:), Bphi_mat(:,:)
        real(kind=8), allocatable :: avg_vals(:,:), timesteps(:)
        type(event), allocatable  :: events_list(:)
        
        integer :: ierr, k, s, idx
        real(kind=8) :: end_time, psi_end
        logical :: flux_av
        integer :: dims_RZ(2)

        ! Only Rank 0 performs this diagnostic logic
        if (sim%my_id /= 0) return

        print *, ">>> Starting Flux Surface Evolution Diagnostic <<<"

        ! 1. Define Flux Surfaces
        allocate(PsiN_list(n_psi_surfaces))
        allocate(Psi_list(n_psi_surfaces))
        do k = 1, n_psi_surfaces
            PsiN_list(k) = 0.001d0 + (k-1) * (0.999d0 - 0.001d0) / real(n_psi_surfaces-1, 8)
            Psi_list(k)  = ES%psi_axis + PsiN_list(k) * (ES%psi_bnd - ES%psi_axis)
        end do

        ! 2. Flux Surface Averages
        allocate(avg_vals(n_psi_surfaces, 20)) ! Adjust size based on JOREK needs
        avg_vals = 0.d0
        flux_av = .true.
        call avg_fluxsurf_list(ES, sim%fields%node_list, sim%fields%element_list, &
                               PsiN_list, avg_vals, flux_av, ierr)
        call save_avg_quantities_h5("flux_averages.h5", PsiN_list, avg_vals)

        ! 3. Geometry and Field Calculation
        call fluxsurface(ES, sim%fields%node_list, sim%fields%element_list, &
                         PsiN_list, R_mat, Z_mat, theta_list, ierr)
        
        dims_RZ = [size(R_mat, 1), size(R_mat, 2)]

        call comp_B_field(ES, sim%fields%node_list, sim%fields%element_list, &
                          R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, ierr)

        ! 4. Save radilal grid for DREAM
        call save_flux_data_h5("flux_surfaces.h5", Psi_list, theta_list, &
                               R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, &
                               ES%R_axis, ES%Z_axis, ES%LCFS_a)

        ! 5. Particle Tracing (Evolution Mapping)
        allocate(particle_fieldline :: sim%groups(1)%particles(dims_RZ(1) * dims_RZ(2)))
        call initialise_on_flux_surfaces(sim, R_mat, Z_mat)

        allocate(events_list(0))
        allocate(timesteps(1))
        timesteps = [1d-6] 
        end_time = sim%time + dt_trace

        call run_particle_trace(sim, timesteps, events_list, end_time)

        ! 6. Mapping Check & Radial Mapping Output
        select type (p => sim%groups(1)%particles)
        type is (particle_fieldline)
            print *, "--- Mapping: psi(n) -> psi(n+1) ---"
            do s = 1, dims_RZ(2)
                idx = s
                call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                    p(idx)%x, psi_end, ierr)
                print "(A,I3,A,F12.6,A,F12.6)", "Surface ", s, " | Init Psi: ", &
                       Psi_list(s), " | End Psi: ", psi_end
            end do
        end select

        ! 7. Clean up
        deallocate(sim%groups(1)%particles, timesteps, events_list)
        deallocate(Psi_list, PsiN_list, theta_list, R_mat, Z_mat)
        deallocate(Br_mat, Bz_mat, Bphi_mat, avg_vals)

        print *, ">>> Flux Surface Evolution Diagnostic Complete <<<"

    end subroutine com_dream

end module mod_com_dream
