
program fifty_steps

    use particle_tracer
    use mod_particle_diagnostics
    use mpi
    use mod_interp
    use mod_atomic_elements
    use mod_particle_evolution
    use mod_particle_recomb
    use mod_particle_conservation
    use mod_particle_io
    use mod_particle_allocation
    use mod_event
    use mod_project_particles
    use mod_jorek_timestepping
    use mod_random_seed
    use mod_basisfunctions
    use nodes_elements
    use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
    use mod_particle_wall_interaction
    use mod_neutral_collision, only: neutral_collisions_from_config, type_neutral_collision
    use mod_projection_functions, only: proj_f_combined_density, proj_f_combined_energy, proj_f_combined_par_momentum
    use mod_particle_puffing
    use mod_edge_domain
    use mod_edge_elements, only: edge_elements
    use mod_atomic_coeff_deuterium, only: ad_deuterium 
    use data_structure, only: type_bnd_element_list, type_bnd_node_list 
    use mod_boundary,   only: boundary_from_grid
    use mod_coupling_settings, only: use_kin_recomb_global
    use mod_initialise_particles
    use equil_info
    use mod_output_file_routines, only: write_to_outputfile
    use mod_expression
    use mod_com_dream, only: com_dream

    use phys_module, only: index_now
    use phys_module, only: tstep,tstep_n,restart_particles, restart, t_start, nout
    use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
    use phys_module, only: n_part_groups, n_aux_var, n_valves_max
    use phys_module, only: nstep_particles, nsubstep_particles, tstep_particles, nout_particles
    use phys_module, only: deuterium_adas,sqrt_mu0_over_rho0
    use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
    use phys_module, only: apply_dirichlet_proj, part_group_configs, init_particles_only
    use phys_module, only: use_manual_random_seed, manual_seed

    use mod_particle_group_id, only: matching_part_config_indices

    use mod_pcg32_rng, only: pcg32_rng
    use mod_rng, only: type_rng, setup_shared_rngs
    implicit none

    ! --- Simulation variables ---
    real(kind=8) :: start_time, end_time, delta_t
    integer :: i_step, inode, k
    integer, parameter :: total_substeps = 50
    real(kind=8)    :: target_time, df
    integer(kind=4) :: ierr
    !logical         :: restart
    type(projection), target                          :: jorek_feedback
    real(kind=8), allocatable :: PsiN_list(:)

    ! Initialize simulation structure
    call sim%initialize()

    if (.true.) then

        sim%time = 23680.550000000017*6.12413767787214d-7

        events = [event(read_jorek_fields_interp_linear(basename='jorek',i=last_file_before_time(sim%time))), & 
                  event(stop_action(), start=26680.550000000017*6.12413767787214d-7)]
                  

        ! Execute first event to load fields into memory at the start time
        call with(sim, events, at=0.0d0)
        !call with(sim, events, at=100.0d0*6.482912854348901d-7)
    endif

    write(*,*) 'Fields loaded. Starting time advancement loop...'

    ! --- Main Event Loop ---

    sim%istep_fluid = index_now

    allocate(PsiN_list(5))
    do k = 1, size(PsiN_list)
        PsiN_list(k) = 0.001d0 + (k-1) * (0.999d0 - 0.001d0) / real(size(PsiN_list)-1, 8)
    end do

    do while (.not. sim%stop_now)

        

        start_time = sim%time
        end_time   = next_event_at(sim, events)
        delta_t    = (end_time - start_time) / real(total_substeps, 8)
        
        write(*,*) 'Advancing interval from ', start_time, ' to ', end_time

        ! 2. Inner loop to call com_dream 50 times
        do i_step = 1, 3

            sim%istep_fluid = sim%istep_fluid + 1
            
            ! Update simulation time for interpolation
            sim%time = start_time + real(i_step, 8) * delta_t
            print *, 'Substep ', i_step, ' at time ', sim%time
            
            ! Execute interpolation event to update fields to the new sim%time
            call with(sim, events, at=sim%time)

            select type (f => sim%fields)
            type is (jorek_fields_interp_linear)
                
                ! Calculate how far we are between time_prev and time_now (0.0 to 1.0)
                df = (f%time_now - sim%time) / (f%time_now - f%time_prev)
                
                !$omp parallel do default(shared) private(inode)
                do inode = 1, f%node_list%n_nodes
                    ! The module defines deltas = Value_now - Value_prev
                    ! To get interpolated value: Value_now - (df * deltas)
                    f%node_list%node(inode)%values = f%node_list%node(inode)%values - df * f%node_list%node(inode)%deltas
                end do
                !$omp end parallel do
            end select

            ! Update equilibrium and projections for the current interpolated state
            if (sim%my_id .eq. 0) then
                call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, &
                                       bnd_node_list, bnd_elm_list, .false.)
            endif
            call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

            call update_equil_state(sim%my_id, sim%fields%node_list, &
                                    sim%fields%element_list, bnd_elm_list, xpoint, xcase )

            jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, &
                                filter_parallel_n0=filter_par_n0, filter = filter_perp, &
                                filter_hyper = filter_hyper, filter_parallel=filter_par, &
                                fractional_digits = 9, do_zonal = .false., &
                                calc_integrals=.false., to_vtk=.false., to_h5 = .false., &
                                basename='projections', nsub=2, do_dirichlet=apply_dirichlet_proj)
            
            aux_node_list => jorek_feedback%node_list

            ! 3. Call com_dream at the current interpolated step
            ! Note: dt_trace might need to be adjusted if it should match the delta_t
            call com_dream(sim, ES, dt_trace= delta_t, PsiN_list=PsiN_list)

        end do 

        ! Ensure the outer loop recognizes we have reached the target
        ! if (sim%time >= next_event_at(sim, events)) then
        !      ! This triggers the stop_action or next file read in the events list
        !      call with(sim, events, at=sim%time)
        ! endif
        sim%stop_now = .true.  

    enddo 

    write(*,*) 'Loop finished at time: ', sim%time

    ! Finalize and clean up
    call sim%finalize

end program fifty_steps
