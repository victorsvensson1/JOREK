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
                               dream_psin, dream_jre_par, dream_bmin
    use libmuscle
    use ymmsl
    use mpi
    implicit none

    private
    public :: com_dream, init_dream_coupling, dream_instance

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


    subroutine com_dream(sim, ES, PsiN_list)
        type(particle_sim), intent(inout) :: sim
        type(t_equil_state), intent(in)   :: ES
        real(kind=8), intent(inout)       :: PsiN_list(:)

        real(kind=8), allocatable :: Psi_list(:), theta_list(:)
        real(kind=8), allocatable :: R_mat(:,:), Z_mat(:,:)
        real(kind=8), allocatable :: Br_mat(:,:), Bz_mat(:,:), Bphi_mat(:,:)
        real(kind=8), allocatable :: avg_vals(:,:), timesteps(:)
        real(kind=8), allocatable :: r_mid_n(:)
        real(kind=8), allocatable :: R_flat(:), Z_flat(:)
        real(kind=8), allocatable :: Br_flat(:), Bz_flat(:), Bphi_flat(:)
        real(kind=8), allocatable :: recv_jre(:), recv_bmin(:)
        type(event), allocatable  :: events_list(:)

        type(LIBMUSCLE_Message)      :: send_msg, recv_msg
        type(LIBMUSCLE_Data)         :: send_data
        type(LIBMUSCLE_Data)         :: d_theta, d_psin, d_T, d_ne, d_Epar, d_rmid
        type(LIBMUSCLE_Data)         :: d_Rax, d_Zax, d_amin, d_dt, d_istep
        type(LIBMUSCLE_Data)         :: d_ntheta, d_Rmat, d_Zmat
        type(LIBMUSCLE_Data)         :: d_Brmat, d_Bzmat, d_Bphimat, d_Psilist
        type(LIBMUSCLE_DataConstRef) :: recv_ref, ref_jre, ref_bmin

        integer :: ierr_m3, ierr_mpi, k, s, n_psi_surfaces, n_theta, n_flat
        real(kind=8) :: psi_end
        logical :: flux_av

        n_psi_surfaces = size(PsiN_list)

        if (sim%my_id /= 0) goto 999

        print *, ">>> com_dream: computing flux surface data <<<"

        allocate(Psi_list(n_psi_surfaces))
        do k = 1, n_psi_surfaces
            Psi_list(k) = ES%psi_axis + PsiN_list(k) * (ES%psi_bnd - ES%psi_axis)
        end do

        flux_av = .true.
        call avg_fluxsurf_list(ES, sim%fields%node_list, sim%fields%element_list, &
                               PsiN_list, avg_vals, flux_av, ierr_m3)
        call fluxsurface(ES, sim%fields%node_list, sim%fields%element_list, &
                         PsiN_list, R_mat, Z_mat, theta_list, ierr_m3)
        call comp_B_field(ES, sim%fields%node_list, sim%fields%element_list, &
                          R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, ierr_m3)

        allocate(r_mid_n(n_psi_surfaces))
        do s = 1, n_psi_surfaces
            r_mid_n(s) = R_mat(1, s) - ES%R_axis
        end do


        if (allocated(sim%groups(1)%particles)) deallocate(sim%groups(1)%particles)
        allocate(particle_gc_relativistic :: sim%groups(1)%particles(size(PsiN_list)))

        call initialise_on_flux_surfaces(sim, R_mat, Z_mat)

        allocate(events_list(0))
        allocate(timesteps(1))
        timesteps = [1d-8]
        call run_particle_trace(sim, timesteps, events_list, sim%time + sim%tstep_fluid_si)

        select type (p => sim%groups(1)%particles)
        type is (particle_gc_relativistic)
            do s = 1, n_psi_surfaces
                call get_psi_at_pos(ES, sim%fields%node_list, sim%fields%element_list, &
                                    p(s)%x, psi_end, ierr_m3)
                PsiN_list(s) = psi_end
            end do
        end select

        n_theta = size(R_mat, 1)
        n_flat  = size(R_mat)
        allocate(R_flat(n_flat), Z_flat(n_flat), Br_flat(n_flat), &
                 Bz_flat(n_flat), Bphi_flat(n_flat))
        R_flat    = reshape(R_mat,    [n_flat])
        Z_flat    = reshape(Z_mat,    [n_flat])
        Br_flat   = reshape(Br_mat,   [n_flat])
        Bz_flat   = reshape(Bz_mat,   [n_flat])
        Bphi_flat = reshape(Bphi_mat, [n_flat])

        d_theta   = LIBMUSCLE_Data_create_grid_1_real8_a(theta_list)
        d_psin    = LIBMUSCLE_Data_create_grid_1_real8_a(PsiN_list)
        d_T       = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,1))
        d_ne      = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,2))
        d_Epar    = LIBMUSCLE_Data_create_grid_1_real8_a(avg_vals(:,3))
        d_rmid    = LIBMUSCLE_Data_create_grid_1_real8_a(r_mid_n)
        d_Rmat    = LIBMUSCLE_Data_create_grid_1_real8_a(R_flat)
        d_Zmat    = LIBMUSCLE_Data_create_grid_1_real8_a(Z_flat)
        d_Brmat   = LIBMUSCLE_Data_create_grid_1_real8_a(Br_flat)
        d_Bzmat   = LIBMUSCLE_Data_create_grid_1_real8_a(Bz_flat)
        d_Bphimat = LIBMUSCLE_Data_create_grid_1_real8_a(Bphi_flat)
        d_Psilist = LIBMUSCLE_Data_create_grid_1_real8_a(Psi_list)
        d_Rax    = LIBMUSCLE_Data_create_real8(ES%R_axis)
        d_Zax    = LIBMUSCLE_Data_create_real8(ES%Z_axis)
        d_amin   = LIBMUSCLE_Data_create_real8(abs(ES%LCFS_a))
        d_dt     = LIBMUSCLE_Data_create_real8(sim%tstep_fluid_si)
        d_istep  = LIBMUSCLE_Data_create_int8(int(sim%istep_fluid, kind=8))
        d_ntheta = LIBMUSCLE_Data_create_int8(int(n_theta, kind=8))

        send_data = LIBMUSCLE_Data_create_dict()
        call LIBMUSCLE_Data_set_item(send_data, 'theta',    d_theta)
        call LIBMUSCLE_Data_set_item(send_data, 'psi_n',    d_psin)
        call LIBMUSCLE_Data_set_item(send_data, 'T',        d_T)
        call LIBMUSCLE_Data_set_item(send_data, 'ne',       d_ne)
        call LIBMUSCLE_Data_set_item(send_data, 'E_par',    d_Epar)
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

        call set_dream_output(psin_arr  = PsiN_list, &
                              jre_arr   = recv_jre,  &
                              bmin_arr  = recv_bmin, &
                              n_pts     = n_dream_psin)

        deallocate(recv_jre, recv_bmin)

        print *, ">>> com_dream: received j_re, n_pts=", n_dream_psin

        deallocate(sim%groups(1)%particles, timesteps, events_list)
        deallocate(Psi_list, theta_list, R_mat, Z_mat)
        deallocate(Br_mat, Bz_mat, Bphi_mat, avg_vals, r_mid_n)
        deallocate(R_flat, Z_flat, Br_flat, Bz_flat, Bphi_flat)

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