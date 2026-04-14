
program victor_test

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
    real(kind=8)    :: target_time
    integer(kind=4) :: ierr
    !logical         :: restart
    type(projection), target                          :: jorek_feedback

    ! Initialize simulation structure
    !call sim%initialize(num_groups=1)
    call sim%initialize()

    !restart = .false. 

    if (.true.) then
        ! Start time
        sim%time = 248680.55000000002 * 6.12413767787214d-7

        ! Initialize events: 
        ! 1. Read fields (linear interpolation)
        ! 2. stop_action at a specific time
        events = [event(read_jorek_fields_interp_linear(basename='jorek',i=last_file_before_time(sim%time))), & 
                  event(stop_action(), start=290680.55000000005*6.12413767787214d-7)]

        ! Execute first event to load fields into memory at the start time
        call with(sim, events, at=0.0d0)
        !call with(sim, events, at=100.0d0*6.482912854348901d-7)
    endif

    write(*,*) 'Fields loaded. Starting time advancement loop...'

    ! --- Main Event Loop ---

    sim%istep_fluid = index_now

    do while (.not. sim%stop_now)

        sim%istep_fluid = index_now
        
        ! Determine the time of the next scheduled event
        target_time = next_event_at(sim, events)
        
        ! Advance simulation time to the target
        sim%time = target_time
        
        write(*,*) 'Advancing to time: ', sim%time

        ! Execute events scheduled for this time (e.g., reading new field files)
        call with(sim, events, at=sim%time)

        if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
        call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

        call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )

        jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0=filter_par_n0,            &
                                filter = filter_perp, filter_hyper = filter_hyper, filter_parallel=filter_par, fractional_digits = 9,       &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections', nsub=2, &
                                do_dirichlet=apply_dirichlet_proj)
        aux_node_list => jorek_feedback%node_list
        

        call com_dream(sim, ES, n_psi_surfaces=50, dt_trace=1.0d-4)



    enddo 

    write(*,*) 'Loop finished at time: ', sim%time

    ! Finalize and clean up
    call sim%finalize

end program victor_test
