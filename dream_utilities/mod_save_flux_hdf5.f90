module mod_save_flux_hdf5
    use hdf5
    use hdf5_io_module
    implicit none

    private 
    public :: save_flux_data_h5, save_avg_quantities_h5

contains

    subroutine save_flux_data_h5(filename, Psi, theta, R, Z, Br, Bz, Bphi, R_ax, Z_ax, a_minor)
        character(len=*), intent(in) :: filename
        real*8, intent(in) :: Psi(:), theta(:), R(:,:), Z(:,:), Br(:,:), Bz(:,:), Bphi(:,:), R_ax, Z_ax, a_minor

#ifdef USE_HDF5
        integer(HID_T) :: i_file, i_group
        integer :: ierr, n_pts, n_psi
        real*8, allocatable :: ptx(:,:), pty(:,:), psi_apRp(:)
        real*8 :: aspect_scale

        ierr = 0
        n_pts = size(R, 1)
        n_psi = size(R, 2)

        allocate(ptx(n_pts, n_psi), pty(n_pts, n_psi), psi_apRp(n_psi))

        ptx = R - R_ax
        pty = Z - Z_ax

        ! Poloidal flux scaling: psi / (R0/a)
        aspect_scale = R_ax / a_minor
        psi_apRp = Psi / aspect_scale

        ! 1. Create the HDF5 file
        call HDF5_create(trim(filename), i_file, ierr)
        if (ierr /= 0) then
            write(*,*) 'Error: Could not create HDF5 file ', trim(filename)
            return
        end if

        ! 2. Create the "equil" group
        call HDF5_group_create(i_file, 'equil'//char(0), i_group)

        ! 3. Save Scalars under group
        call HDF5_real_saving(i_group, R_ax, 'Rp'//char(0))
        call HDF5_real_saving(i_group, Z_ax, 'Zp'//char(0))

        ! 4. Save 1D Arrays under group
        call HDF5_array1D_saving(i_group, psi_apRp, n_psi, 'psi_apRp'//char(0))
        call HDF5_array1D_saving(i_group, theta, n_pts, 'theta'//char(0))

        ! 5. Save 2D Matrices under group
        call HDF5_array2D_saving(i_group, transpose(ptx), n_psi, n_pts, 'ptx'//char(0))
        call HDF5_array2D_saving(i_group, transpose(pty),    n_psi, n_pts, 'pty'//char(0))
        call HDF5_array2D_saving(i_group, transpose(Br),   n_psi, n_pts, 'ptBx'//char(0))
        call HDF5_array2D_saving(i_group, transpose(Bz),   n_psi, n_pts, 'ptBy'//char(0))
        call HDF5_array2D_saving(i_group, transpose(Bphi), n_psi, n_pts, 'ptBPHI'//char(0))

        ! 6. Close group and file
        call HDF5_group_close(i_group)
        call HDF5_close(i_file)
        deallocate(ptx, pty, psi_apRp)
#else
        write(*,*) 'Warning: HDF5 not enabled. File not saved.'
#endif
    end subroutine save_flux_data_h5


    subroutine save_avg_quantities_h5(filename, Psi, avg_vals)
        character(len=*), intent(in) :: filename
        real*8, intent(in) :: Psi(:)
        real*8, intent(in) :: avg_vals(:,:)

#ifdef USE_HDF5
        integer(HID_T) :: i_file, i_group
        integer :: ierr, n_psi

        ierr = 0
        n_psi = size(Psi)

        call HDF5_create(trim(filename), i_file, ierr)
        if (ierr /= 0) return

        ! Use group for average quantities as well
        call HDF5_group_create(i_file, 'equil'//char(0), i_group)

        call HDF5_array1D_saving(i_group, Psi, n_psi, 'Psi'//char(0))
        call HDF5_array1D_saving(i_group, avg_vals(:,1), n_psi, 'T'//char(0))
        call HDF5_array1D_saving(i_group, avg_vals(:,2), n_psi, 'ne'//char(0))
        call HDF5_array1D_saving(i_group, avg_vals(:,3), n_psi, 'E_parallel'//char(0))

        call HDF5_group_close(i_group)
        call HDF5_close(i_file)
#else
        write(*,*) 'Warning: HDF5 not enabled. Average data not saved.'
#endif
    end subroutine save_avg_quantities_h5

    ! -----------------------------------------------------------
    ! Internal Private Helper Subroutines
    ! -----------------------------------------------------------

    subroutine HDF5_group_create(file_id, group_name, group_id)
        integer(HID_T), intent(in)  :: file_id
        character(len=*), intent(in):: group_name
        integer(HID_T), intent(out) :: group_id
        integer :: ierr
        ! Standard HDF5 Fortran API call
        call h5gcreate_f(file_id, trim(group_name), group_id, ierr)
    end subroutine HDF5_group_create

    subroutine HDF5_group_close(group_id)
        integer(HID_T), intent(in) :: group_id
        integer :: ierr
        ! Standard HDF5 Fortran API call
        call h5gclose_f(group_id, ierr)
    end subroutine HDF5_group_close

end module mod_save_flux_hdf5