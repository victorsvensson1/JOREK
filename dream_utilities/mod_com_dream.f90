module mod_com_dream
    use particle_tracer
    use mod_fluxsurf_evol
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

    subroutine com_dream(sim, ES, dt_trace, PsiN_list)
        type(particle_sim), intent(inout) :: sim
        type(t_equil_state), intent(in)   :: ES
        real(kind=8), intent(in)          :: dt_trace
        real(kind=8), intent(inout)          :: PsiN_list(:)

        ! --- Local Data ---
        real(kind=8), allocatable :: Psi_list(:), theta_list(:)
        real(kind=8), allocatable :: R_mat(:,:), Z_mat(:,:)
        real(kind=8), allocatable :: Br_mat(:,:), Bz_mat(:,:), Bphi_mat(:,:)
        real(kind=8), allocatable :: avg_vals(:,:), timesteps(:)
        type(event), allocatable  :: events_list(:)
        real(kind=8), allocatable :: r_mid_n(:)
        
        integer :: ierr, k, s, idx, n_psi_surfaces
        real(kind=8) :: end_time, psi_end
        logical :: flux_av
        character(len=64) :: fname_avg, fname_surf, fname_map
        real(kind=8) :: R_new

        n_psi_surfaces = size(PsiN_list)

        ! Only Rank 0 performs this diagnostic logic
        if (sim%my_id /= 0) return

        print *, ">>> Starting Flux Surface Evolution Diagnostic <<<"

        ! 1. Define Flux Surfaces
        allocate(Psi_list(n_psi_surfaces))
        
        do k = 1, n_psi_surfaces
            !PsiN_list(k) = 0.001d0 + (k-1) * (0.999d0 - 0.001d0) / real(n_surfaces-1, 8)
            Psi_list(k)  = ES%psi_axis + PsiN_list(k) * (ES%psi_bnd - ES%psi_axis)
        end do

        print *, "PsiN_list (normalized): ", PsiN_list

        ! 2. Flux Surface Averages
        flux_av = .true.
        ierr = 0
        call avg_fluxsurf_list(ES, sim%fields%node_list, sim%fields%element_list, &
                               PsiN_list, avg_vals, flux_av, ierr)
        !call avg_fluxsurf_list(avg_vals, n_psi_surfaces)
        write(fname_avg, '("flux_averages_", I6.6, ".h5")') sim%istep_fluid
        call save_avg_quantities_h5(trim(fname_avg), PsiN_list, avg_vals)

        ! 3. Geometry and Field Calculation
        call fluxsurface(ES, sim%fields%node_list, sim%fields%element_list, &
                         PsiN_list, R_mat, Z_mat, theta_list, ierr)

        call comp_B_field(ES, sim%fields%node_list, sim%fields%element_list, &
                          R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, ierr)

        ! 4. Save radilal grid for DREAM
        write(fname_surf, '("flux_surfaces_", I6.6, ".h5")') sim%istep_fluid
        call save_flux_data_h5(trim(fname_surf), Psi_list, theta_list, &
                               R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, &
                               ES%R_axis, ES%Z_axis, ES%LCFS_a)

        ! 5. Particle Tracing (Evolution Mapping)
        allocate(r_mid_n(n_psi_surfaces))
        do s = 1, n_psi_surfaces
            r_mid_n(s) = R_mat(1, s) - ES%R_axis
        end do

        allocate(particle_gc_relativistic :: sim%groups(1)%particles(size(PsiN_list)))
        sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU
        call initialise_on_flux_surfaces(sim, R_mat, Z_mat)
        do s = 1, n_psi_surfaces
            print *, "Initialized particle ", s, " at R=", sim%groups(1)%particles(s)%x(1), &
                     " Z=", sim%groups(1)%particles(s)%x(2)
        end do

        allocate(events_list(0))
        allocate(timesteps(1))
        timesteps = [1d-8] 
        end_time = sim%time + dt_trace

        call run_particle_trace(sim, timesteps, events_list, end_time)

        ! 6. Mapping Check & Radial Mapping Output
        select type (p => sim%groups(1)%particles)
        type is (particle_gc_relativistic)
            print *, "--- Mapping: psi(n) -> psi(n+1) ---"
            do s = 1, n_psi_surfaces
                call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                    p(s)%x, psi_end, ierr)
                print "(A,I3,A,F12.6,A,F12.6)", "Surface ", s, " | Init Psi: ", &
                       PsiN_list(s), " | End Psi: ", psi_end
                
                PsiN_list(s) = psi_end
            end do
        end select
        
        write(fname_map, '("radial_mapping_", I6.6, ".h5")') sim%istep_fluid
        call save_radial_mapping_h5(trim(fname_map), r_mid_n)

        ! 7. Clean up
        deallocate(sim%groups(1)%particles, timesteps, events_list)
        deallocate(Psi_list, theta_list, R_mat, Z_mat)
        deallocate(Br_mat, Bz_mat, Bphi_mat, avg_vals)
        deallocate(r_mid_n)

        print *, ">>> Flux Surface Evolution Diagnostic Complete <<<"

    end subroutine com_dream


end module mod_com_dream
