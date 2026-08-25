module mod_com_dream
    use particle_tracer
    use mod_fluxsurf_evol
    use mod_computeB
    use mod_fields_linear
    use equil_info
    use mod_fluxsurf_evol, only: initialise_on_flux_surfaces, run_particle_trace, get_psi_at_pos
    use mod_fluxsurf_avg, only: avg_fluxsurf_list
    use mod_fluxsurf_compute, only: fluxsurface
    use mod_computeB, only: comp_B_field
    use mod_save_flux_hdf5
    use mod_dream_input, only: set_dream_output, n_dream_psin, &
                               dream_psin, dream_jre_par, dream_bmin, &
                               load_dream_output_from_file
    use libmuscle
    use ymmsl
    use mpi
    implicit none

    private
    public :: com_dream, init_dream_coupling, restart_dream_coupling, dream_instance

    type(LIBMUSCLE_Instance), save :: dream_instance
    logical, save :: dream_particles_allocated = .false.

contains

    subroutine init_dream_coupling()
        type(LIBMUSCLE_PortsDescription) :: ports
        integer :: my_id, ierr

        call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
        if (my_id == 0) then
            ports = LIBMUSCLE_PortsDescription()
            call ports%add(YMMSL_Operator_O_I, 'plasma_state_out')
            call ports%add(YMMSL_Operator_S,   'jre_in')
            dream_instance = LIBMUSCLE_Instance(ports)
            call LIBMUSCLE_PortsDescription_free(ports)
        end if
    end subroutine init_dream_coupling

     subroutine restart_dream_coupling()

        character(len=1024) :: restart_file
        integer :: my_id, ierr, ierr_mpi, len_env, stat_env

        call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

        if (my_id == 0) then
            call get_environment_variable('JOREK_DREAM_RESTART_FILE', restart_file, len_env, stat_env)
            if (stat_env == 0 .and. len_env > 0) then
                write(*,*) '>>> restart_dream_coupling: loading ', trim(restart_file)
                call load_dream_output_from_file(trim(restart_file), ierr)
                if (ierr /= 0) then
                    write(*,*) 'ERROR: restart_dream_coupling failed to load ', trim(restart_file)
                    n_dream_psin = 0
                else
                    write(*,*) '>>> restart_dream_coupling: primed with n_pts=', n_dream_psin
                end if
            else
                write(*,*) '>>> restart_dream_coupling: JOREK_DREAM_RESTART_FILE not set (stat=', &
                           stat_env, ', len=', len_env, ') — starting with empty DREAM cache'
                n_dream_psin = 0
            end if
        end if

        call MPI_Bcast(n_dream_psin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
        if (n_dream_psin > 0) then
            if (my_id /= 0) then
                if (allocated(dream_psin))    deallocate(dream_psin)
                if (allocated(dream_jre_par)) deallocate(dream_jre_par)
                if (allocated(dream_bmin))    deallocate(dream_bmin)
                allocate(dream_psin(n_dream_psin))
                allocate(dream_jre_par(n_dream_psin))
                allocate(dream_bmin(n_dream_psin))
            end if
            call MPI_Bcast(dream_psin,    n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
            call MPI_Bcast(dream_jre_par, n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
            call MPI_Bcast(dream_bmin,    n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
        end if
    end subroutine restart_dream_coupling


    subroutine com_dream(sim, ES, n_psi_surfaces)
        type(particle_sim), intent(inout) :: sim
        type(t_equil_state), intent(in)   :: ES
        integer, intent(in)               :: n_psi_surfaces

        real(kind=8), allocatable :: avg_vals(:,:), timesteps(:)
        real(kind=8), allocatable :: r_mid_n(:)
        real(kind=8), allocatable :: recv_jre(:), recv_bmin(:)
        type(event), allocatable  :: events_list(:)

        ! --- Fixed uniform psi_n grid, used for EVERYTHING sent to DREAM:
        ! the equilibrium geometry (LUKE file) AND the flux-surface-averaged
        ! T/ne/E_par, so the two stay consistent with each other and DREAM's
        ! domain always spans the full closed-flux-surface region.
        real(kind=8), allocatable :: PsiN_geom(:), Psi_list_geom(:), theta_list_geom(:)
        real(kind=8), allocatable :: R_mat_geom(:,:), Z_mat_geom(:,:)
        real(kind=8), allocatable :: Br_mat_geom(:,:), Bz_mat_geom(:,:), Bphi_mat_geom(:,:)
        real(kind=8), allocatable :: R_flat_geom(:), Z_flat_geom(:)
        real(kind=8), allocatable :: Br_flat_geom(:), Bz_flat_geom(:), Bphi_flat_geom(:)
        integer :: n_theta_geom, n_flat_geom, ierr_geom

        real(kind=8), allocatable :: psin_mapped(:)
        type(LIBMUSCLE_Data)      :: d_psin_mapped

        type(LIBMUSCLE_Message)      :: send_msg, recv_msg
        type(LIBMUSCLE_Data)         :: send_data
        type(LIBMUSCLE_Data)         :: d_theta, d_psin, d_T, d_ne, d_Epar, d_rmid, d_Zeff
        type(LIBMUSCLE_Data)         :: d_Rax, d_Zax, d_amin, d_dt, d_istep
        type(LIBMUSCLE_Data)         :: d_ntheta, d_Rmat, d_Zmat
        type(LIBMUSCLE_Data)         :: d_Brmat, d_Bzmat, d_Bphimat, d_Psilist
        type(LIBMUSCLE_DataConstRef) :: recv_ref, ref_jre, ref_bmin

        integer :: ierr_m3, ierr_mpi, k, s
        real(kind=8) :: psi_end
        logical :: flux_av

        type(particle_gc_relativistic), allocatable :: trace_particles(:)
        real(kind=8), parameter :: electron_mass = 5.4857990907016d-4

        logical           :: disable_tracing
        character(len=8)  :: env_disable_tracing
        integer           :: len_env_dt, stat_env_dt
        
        !for conv test
        logical           :: do_convergence_check
        character(len=8)  :: env_conv_check
        character(len=16) :: env_conv_levels
        integer           :: len_env_cc, stat_env_cc, len_env_cl, stat_env_cl
        integer           :: n_conv_levels, lvl
        real(kind=8)      :: x_factor, t_end
        type(particle_gc_relativistic), allocatable :: trace_particles_conv(:)
        real(kind=8), allocatable :: psin_mapped_conv(:)


        logical           :: do_dt_convergence_check
        character(len=8)  :: env_dt_check
        integer           :: len_env_dtc, stat_env_dtc
        integer           :: idx
        real(kind=8), parameter :: dt_levels(9) = &
            [1d-7, 5d-8, 2.5d-8, 1.25d-8, 6.25d-9, 3.125d-9, 1.5625d-9, 7.8125d-10, 3.90625d-10]
        type(particle_gc_relativistic), allocatable :: trace_particles_dt(:)
        real(kind=8), allocatable :: psin_mapped_dt(:)
        real(kind=8), allocatable :: dt_step(:)
        ! end conv test

        if (sim%my_id /= 0) goto 999

        call get_environment_variable('JOREK_DREAM_DISABLE_TRACING', env_disable_tracing, len_env_dt, stat_env_dt)
        disable_tracing = (stat_env_dt == 0 .and. len_env_dt > 0 .and. trim(env_disable_tracing) == '1')

        ! conv tracing
        call get_environment_variable('JOREK_DREAM_DISABLE_TRACING', env_disable_tracing, len_env_dt, stat_env_dt)
        disable_tracing = (stat_env_dt == 0 .and. len_env_dt > 0 .and. trim(env_disable_tracing) == '1')

        ! --- Convergence-check flag: re-traces with the interval shrunk by
        ! successive factors of 2 and prints the mapped psi_N at each level,
        ! purely as a diagnostic. Does not affect the value sent to DREAM.
        call get_environment_variable('JOREK_DREAM_CONVERGENCE_CHECK', env_conv_check, len_env_cc, stat_env_cc)
        do_convergence_check = (stat_env_cc == 0 .and. len_env_cc > 0 .and. trim(env_conv_check) == '1')

        n_conv_levels = 5   ! default: x = 1, 1/2, 1/4, 1/8, 1/16
        call get_environment_variable('JOREK_DREAM_CONVERGENCE_LEVELS', env_conv_levels, len_env_cl, stat_env_cl)
        if (stat_env_cl == 0 .and. len_env_cl > 0) read(env_conv_levels, *) n_conv_levels

        call get_environment_variable('JOREK_DREAM_TIMESTEP_CHECK', env_dt_check, len_env_dtc, stat_env_dtc)
        do_dt_convergence_check = (stat_env_dtc == 0 .and. len_env_dtc > 0 .and. trim(env_dt_check) == '1')

        ! conv tracing end

        print *, ">>> com_dream: computing flux surface data <<<"

        allocate(PsiN_geom(n_psi_surfaces))
        do k = 1, n_psi_surfaces
            PsiN_geom(k) = 0.001d0 + (k-1) * (0.999d0 - 0.001d0) / real(n_psi_surfaces-1, 8)
        end do

        allocate(Psi_list_geom(n_psi_surfaces))
        do k = 1, n_psi_surfaces
            Psi_list_geom(k) = ES%psi_axis + PsiN_geom(k) * (ES%psi_bnd - ES%psi_axis)
        end do

        flux_av = .true.
        call avg_fluxsurf_list(ES, sim%fields%node_list, sim%fields%element_list, &
                               PsiN_geom, avg_vals, flux_av, ierr_m3)
        call fluxsurface(ES, sim%fields%node_list, sim%fields%element_list, &
                         PsiN_geom, R_mat_geom, Z_mat_geom, theta_list_geom, ierr_geom)
        call comp_B_field(ES, sim%fields%node_list, sim%fields%element_list, &
                          R_mat_geom, Z_mat_geom, Br_mat_geom, Bz_mat_geom, Bphi_mat_geom, ierr_geom)

        allocate(r_mid_n(n_psi_surfaces))
        do s = 1, n_psi_surfaces
            r_mid_n(s) = R_mat_geom(1, s) - ES%R_axis
        end do

        ! Administrative allocations kept unconditional (cheap either way, and
        ! needed regardless by the deallocate() call at the end of this routine).
        allocate(trace_particles(n_psi_surfaces))
        allocate(events_list(0))
        allocate(timesteps(1))
        timesteps = [1d-9]

        if (disable_tracing) then

            allocate(psin_mapped(n_psi_surfaces))
            psin_mapped = PsiN_geom
            print *, ">>> com_dream: flux-surface tracing DISABLED ", &
                     "(JOREK_DREAM_DISABLE_TRACING=1) -- psin_mapped = PsiN_geom (identity) <<<"

        else

            ! Seed markers on the FIXED grid (same geometry just built above),
            ! push them across this fluid step, and read off where each ends up
            ! in normalized flux -- this is psin_mapped, the per-step mapping.
            call initialise_on_flux_surfaces(sim, trace_particles, R_mat_geom, Z_mat_geom)

            if (sim%istep_fluid > 1) then
                select type (fields => sim%fields)
                type is (jorek_fields_interp_linear)
                    call run_particle_trace(sim, trace_particles, electron_mass, timesteps, &
                                            fields%time_prev, fields%time_now)
                end select
            end if

            allocate(psin_mapped(n_psi_surfaces))
            do s = 1, n_psi_surfaces
                call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                    trace_particles(s)%x, psi_end, ierr_m3)
                psin_mapped(s) = psi_end
            end do

            print *, ">>> com_dream: flux-surface tracing results (start psiN -> mapped psiN) <<<"
            do s = 1, n_psi_surfaces
                write(*,'(A,I5,A,F18.14,A,F18.14)') '    s=', s, &
                    '  psiN_start=', PsiN_geom(s), &
                    '  psiN_mapped=', psin_mapped(s)
            end do

            ! --- Optional convergence check: shrink the tracing interval by
            ! successive factors of 2 and compare against the full-interval
            ! result computed above. Diagnostic only -- psin_mapped (used
            ! below to send data to DREAM) is left untouched.
            if (do_convergence_check .and. sim%istep_fluid > 1) then

                print *, ">>> com_dream: convergence check active " // &
                         "(JOREK_DREAM_CONVERGENCE_CHECK=1), n_levels=", n_conv_levels, " <<<"

                allocate(trace_particles_conv(n_psi_surfaces))
                allocate(psin_mapped_conv(n_psi_surfaces))

                do lvl = 1, n_conv_levels
                    x_factor = 1d0 / 2d0**(lvl-1)

                    select type (fields => sim%fields)
                    type is (jorek_fields_interp_linear)
                        t_end = fields%time_prev + x_factor * (fields%time_now - fields%time_prev)

                        ! re-seed from the same starting flux surfaces every time --
                        ! run_particle_trace advances trace_particles_conv%x in place,
                        ! so each level must start fresh from the unperturbed geometry.
                        call initialise_on_flux_surfaces(sim, trace_particles_conv, R_mat_geom, Z_mat_geom)

                        call run_particle_trace(sim, trace_particles_conv, electron_mass, timesteps, &
                                                fields%time_prev, t_end)
                    end select

                    do s = 1, n_psi_surfaces
                        call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                            trace_particles_conv(s)%x, psi_end, ierr_m3)
                        psin_mapped_conv(s) = psi_end
                    end do

                    print *, '    --- level ', lvl, ': x=', x_factor, &
                             '  t_end=', t_end, ' ---'
                    do s = 1, n_psi_surfaces
                        write(*,'(A,I5,A,F21.17,A,F21.17,A,ES25.17E3)') '        s=', s, &
                            '  psiN_start=', PsiN_geom(s), &
                            '  psiN_mapped=', psin_mapped_conv(s), &
                            '  diff_vs_full=', psin_mapped_conv(s) - psin_mapped(s)
                    end do
                end do

                deallocate(trace_particles_conv, psin_mapped_conv)

            end if

            ! --- Optional timestep convergence check: trace across the FULL
            ! interval [time_prev, time_now] (x=1, unchanged) but shrink the
            ! internal push step size, to see whether the mapped psi_N has
            ! converged with respect to the integrator's step size itself
            ! (independent of the interval-shrinking check above).
            if (do_dt_convergence_check .and. sim%istep_fluid > 1) then

                print *, ">>> com_dream: timestep convergence check active " // &
                         "(JOREK_DREAM_TIMESTEP_CHECK=1), n_levels=", size(dt_levels), " <<<"

                allocate(trace_particles_dt(n_psi_surfaces))
                allocate(psin_mapped_dt(n_psi_surfaces))
                allocate(dt_step(1))

                do idx = 1, size(dt_levels)
                    dt_step(1) = dt_levels(idx)

                    select type (fields => sim%fields)
                    type is (jorek_fields_interp_linear)
                        ! re-seed from the same starting flux surfaces every
                        ! time, same as the interval-convergence loop
                        call initialise_on_flux_surfaces(sim, trace_particles_dt, R_mat_geom, Z_mat_geom)

                        call run_particle_trace(sim, trace_particles_dt, electron_mass, dt_step, &
                                                fields%time_prev, fields%time_now)
                    end select

                    do s = 1, n_psi_surfaces
                        call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                            trace_particles_dt(s)%x, psi_end, ierr_m3)
                        psin_mapped_dt(s) = psi_end
                    end do

                    print *, '    --- dt level ', idx, ': dt=', dt_step(1), ' ---'
                    do s = 1, n_psi_surfaces
                        write(*,'(A,I5,A,F21.17,A,F21.17,A,ES25.17E3)') '        s=', s, &
                            '  psiN_start=', PsiN_geom(s), &
                            '  psiN_mapped=', psin_mapped_dt(s), &
                            '  diff_vs_full=', psin_mapped_dt(s) - psin_mapped(s)
                    end do
                end do

                deallocate(trace_particles_dt, psin_mapped_dt, dt_step)

            end if

        end if


        n_theta_geom = size(R_mat_geom, 1)
        n_flat_geom  = size(R_mat_geom)
        allocate(R_flat_geom(n_flat_geom), Z_flat_geom(n_flat_geom), Br_flat_geom(n_flat_geom), &
                 Bz_flat_geom(n_flat_geom), Bphi_flat_geom(n_flat_geom))
        R_flat_geom    = reshape(R_mat_geom,    [n_flat_geom])
        Z_flat_geom    = reshape(Z_mat_geom,    [n_flat_geom])
        Br_flat_geom   = reshape(Br_mat_geom,   [n_flat_geom])
        Bz_flat_geom   = reshape(Bz_mat_geom,   [n_flat_geom])
        Bphi_flat_geom = reshape(Bphi_mat_geom, [n_flat_geom])

        d_theta   = LIBMUSCLE_Data_create_grid_1_real8_a(theta_list_geom)
        d_psin    = LIBMUSCLE_Data_create_grid_1_real8_a(PsiN_geom)
        d_T       = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,1))
        d_ne      = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,2))
        d_Epar    = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,3))
        d_Zeff    = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,4))
        d_rmid    = LIBMUSCLE_Data_create_grid_1_real8_a(r_mid_n)
        d_Rmat    = LIBMUSCLE_Data_create_grid_1_real8_a(R_flat_geom)
        d_Zmat    = LIBMUSCLE_Data_create_grid_1_real8_a(Z_flat_geom)
        d_Brmat   = LIBMUSCLE_Data_create_grid_1_real8_a(Br_flat_geom)
        d_Bzmat   = LIBMUSCLE_Data_create_grid_1_real8_a(Bz_flat_geom)
        d_Bphimat = LIBMUSCLE_Data_create_grid_1_real8_a(Bphi_flat_geom)
        d_Psilist = LIBMUSCLE_Data_create_grid_1_real8_a(Psi_list_geom)
        d_psin_mapped = LIBMUSCLE_Data_create_grid_1_real8_a(psin_mapped)
        d_Rax    = LIBMUSCLE_Data_create_real8(ES%R_axis)
        d_Zax    = LIBMUSCLE_Data_create_real8(ES%Z_axis)
        d_amin   = LIBMUSCLE_Data_create_real8(abs(ES%LCFS_a))
        d_dt     = LIBMUSCLE_Data_create_real8(sim%tstep_fluid_si)
        d_istep  = LIBMUSCLE_Data_create_int8(int(sim%istep_fluid, kind=8))
        d_ntheta = LIBMUSCLE_Data_create_int8(int(n_theta_geom, kind=8))

        send_data = LIBMUSCLE_Data_create_dict()
        call LIBMUSCLE_Data_set_item(send_data, 'theta',    d_theta)
        call LIBMUSCLE_Data_set_item(send_data, 'psi_n',    d_psin)
        call LIBMUSCLE_Data_set_item(send_data, 'T',        d_T)
        call LIBMUSCLE_Data_set_item(send_data, 'ne',       d_ne)
        call LIBMUSCLE_Data_set_item(send_data, 'E_par',    d_Epar)
        call LIBMUSCLE_Data_set_item(send_data, 'Zeff',     d_Zeff)
        call LIBMUSCLE_Data_set_item(send_data, 'r_mid_n',  d_rmid)
        call LIBMUSCLE_Data_set_item(send_data, 'R_axis',   d_Rax)
        call LIBMUSCLE_Data_set_item(send_data, 'Z_axis',   d_Zax)
        call LIBMUSCLE_Data_set_item(send_data, 'a_minor',  d_amin)
        call LIBMUSCLE_Data_set_item(send_data, 'dt',       d_dt)
        call LIBMUSCLE_Data_set_item(send_data, 'istep',    d_istep)
        call LIBMUSCLE_Data_set_item(send_data, 'n_theta',  d_ntheta)
        call LIBMUSCLE_Data_set_item(send_data, 'R_mat',    d_Rmat)
        call LIBMUSCLE_Data_set_item(send_data, 'Z_mat',    d_Zmat)
        call LIBMUSCLE_Data_set_item(send_data, 'Br_mat',   d_Brmat)
        call LIBMUSCLE_Data_set_item(send_data, 'Bz_mat',   d_Bzmat)
        call LIBMUSCLE_Data_set_item(send_data, 'Bphi_mat', d_Bphimat)
        call LIBMUSCLE_Data_set_item(send_data, 'Psi_list', d_Psilist)
        call LIBMUSCLE_Data_set_item(send_data, 'psi_n_mapped', d_psin_mapped)

        send_msg = LIBMUSCLE_Message_create(sim%time, send_data)
        call LIBMUSCLE_Instance_send(dream_instance, 'plasma_state_out', send_msg)

        call LIBMUSCLE_Message_free(send_msg)
        call LIBMUSCLE_Data_free(send_data)
        call LIBMUSCLE_Data_free(d_theta);   call LIBMUSCLE_Data_free(d_psin)
        call LIBMUSCLE_Data_free(d_T);       call LIBMUSCLE_Data_free(d_ne)
        call LIBMUSCLE_Data_free(d_Epar);    call LIBMUSCLE_Data_free(d_rmid)
        call LIBMUSCLE_Data_free(d_Rax);     call LIBMUSCLE_Data_free(d_Zax)
        call LIBMUSCLE_Data_free(d_amin);    call LIBMUSCLE_Data_free(d_dt)
        call LIBMUSCLE_Data_free(d_istep);   call LIBMUSCLE_Data_free(d_ntheta)
        call LIBMUSCLE_Data_free(d_Rmat);    call LIBMUSCLE_Data_free(d_Zmat)
        call LIBMUSCLE_Data_free(d_Brmat);   call LIBMUSCLE_Data_free(d_Bzmat)
        call LIBMUSCLE_Data_free(d_Bphimat); call LIBMUSCLE_Data_free(d_Psilist)
        call LIBMUSCLE_Data_free(d_psin_mapped)
        call LIBMUSCLE_Data_free(d_Zeff)

        ! 8. Receive j_re and Bmin from DREAM
        recv_msg = LIBMUSCLE_Instance_receive(dream_instance, 'jre_in')
        recv_ref = LIBMUSCLE_Message_get_data(recv_msg)

        ref_jre  = LIBMUSCLE_DataConstRef_get_item_by_key(recv_ref, 'jre_par')
        ref_bmin = LIBMUSCLE_DataConstRef_get_item_by_key(recv_ref, 'bmin')

        n_dream_psin = int(LIBMUSCLE_DataConstRef_size(ref_jre))

        allocate(recv_jre(n_dream_psin), recv_bmin(n_dream_psin))
        call LIBMUSCLE_DataConstRef_elements(ref_jre,  recv_jre,  ierr_m3)
        call LIBMUSCLE_DataConstRef_elements(ref_bmin, recv_bmin, ierr_m3)

        call LIBMUSCLE_Message_free(recv_msg)

        ! recv_jre/recv_bmin come back defined on DREAM's own radial grid,
        ! which was built from PsiN_geom (sent as 'psi_n' above) -- so the
        ! cache used by interp_dream_jre to place j_re back onto JOREK's mesh
        ! must be keyed by PsiN_geom too.
        call set_dream_output(psin_arr  = PsiN_geom, &
                              jre_arr   = recv_jre,  &
                              bmin_arr  = recv_bmin, &
                              n_pts     = n_dream_psin)

        deallocate(recv_jre, recv_bmin)

        print *, ">>> com_dream: received j_re, n_pts=", n_dream_psin

        deallocate(trace_particles, timesteps, events_list)
        deallocate(avg_vals, r_mid_n, psin_mapped)
        deallocate(PsiN_geom, Psi_list_geom, theta_list_geom, R_mat_geom, Z_mat_geom)
        deallocate(Br_mat_geom, Bz_mat_geom, Bphi_mat_geom)
        deallocate(R_flat_geom, Z_flat_geom, Br_flat_geom, Bz_flat_geom, Bphi_flat_geom)

        print *, ">>> com_dream complete <<<"

999 continue

        call MPI_Bcast(n_dream_psin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
        if (n_dream_psin > 0) then
            if (sim%my_id /= 0) then
                if (allocated(dream_psin))    deallocate(dream_psin)
                if (allocated(dream_jre_par)) deallocate(dream_jre_par)
                if (allocated(dream_bmin))    deallocate(dream_bmin)
                allocate(dream_psin(n_dream_psin))
                allocate(dream_jre_par(n_dream_psin))
                allocate(dream_bmin(n_dream_psin))
            end if
            call MPI_Bcast(dream_psin,    n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
            call MPI_Bcast(dream_jre_par, n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
            call MPI_Bcast(dream_bmin,    n_dream_psin, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
        end if

    end subroutine com_dream

end module mod_com_dream