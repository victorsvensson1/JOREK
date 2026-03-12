!> Module for particle tracing in jorek to keep track of flux surfaces.



module mod_particle_tracing

    use particle_tracer

    use mod_gc_relativistic

    implicit none



contains



    subroutine run_particle_trace(sim, timesteps, events)

        type(particle_sim), intent(inout) :: sim

        real(kind=8), intent(in)                 :: timesteps(:)

        type(event), intent(inout)               :: events(:)

        

        real(kind=8)    :: target_time

        integer(kind=4) :: i, j, k, n_steps, n_lost



        ! Loop until the simulation is stopped (triggered by stop_action event)

        do while (.not. sim%stop_now)

            ! Extract the next event time (e.g., next diagnostic or stop time)

            target_time = next_event_at(sim, events)

            

            do i = 1, size(sim%groups)

                n_steps = nint((target_time - sim%time) / timesteps(i))

                n_lost = 0



                select type (particles => sim%groups(i)%particles)

                type is (particle_gc_relativistic)  

                    do j = 1, size(particles, 1)

                        do k = 1, n_steps

                            if (particles(j)%i_elm <= 0) exit



                            sim%time = sim%time + timesteps(i)

                            

                            ! Perform the actual push in JOREK fields

                            call runge_kutta_fixed_dt_gc_push_jorek( &

                                sim%fields, sim%time, timesteps(i), &

                                sim%groups(i)%mass, particles(j))



                            if (particles(j)%i_elm <= 0) then

                                n_lost = n_lost + 1

                                write(*,*) 'PARTICLE LOST AT TIME: ', sim%time

                                sim%stop_now = .true. 

                                exit

                            end if    

                        end do 

                    end do

                end select

            end do



            ! Update current time and run scheduled events (diagnostics, etc.)

            sim%time = target_time

            call with(sim, events, at=sim%time)

        end do

    end subroutine run_particle_trace



end module mod_particle_tracing
