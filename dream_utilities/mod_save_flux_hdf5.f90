module mod_save_flux_hdf5
    use hdf5
    use hdf5_io_module
    implicit none

    private 
    public :: save_flux_data_h5, save_avg_quantities_h5

contains
    subroutine save_flux_data_h5(filename, Psi, theta, R, Z, Br, Bz, Bphi, R_ax, Z_ax)


    character(len=*), intent(in) :: filename
    real*8, intent(in) :: Psi(:), theta(:), R(:,:), Z(:,:), Br(:,:), Bz(:,:), Bphi(:,:), R_ax, Z_ax

#ifdef USE_HDF5
    integer(HID_T) :: i_file
    integer :: ierr, n_pts, n_psi

    ierr = 0
    n_pts = size(R, 1)
    n_psi = size(R, 2)

    ! 1. Create the HDF5 file using your custom wrapper
    call HDF5_create(trim(filename), i_file, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error: Could not create HDF5 file ', trim(filename)
        return
    end if

    ! 2. Save Scalars (Axis positions)
    call HDF5_real_saving(i_file, R_ax, 'R_axis'//char(0))
    call HDF5_real_saving(i_file, Z_ax, 'Z_axis'//char(0))

    ! 3. Save 1D Arrays (Psi and Theta)
    call HDF5_array1D_saving(i_file, Psi, n_psi, 'Psi_list'//char(0))
    call HDF5_array1D_saving(i_file, theta, n_pts, 'theta_list'//char(0))

    ! 4. Save 2D Matrices (R, Z, and B-fields)
    ! Note: Your wrapper likely expects (file_id, array, dim1, dim2, name)
    call HDF5_array2D_saving(i_file, R,    n_pts, n_psi, 'R_mat'//char(0))
    call HDF5_array2D_saving(i_file, Z,    n_pts, n_psi, 'Z_mat'//char(0))
    call HDF5_array2D_saving(i_file, Br,   n_pts, n_psi, 'Br_mat'//char(0))
    call HDF5_array2D_saving(i_file, Bz,   n_pts, n_psi, 'Bz_mat'//char(0))
    call HDF5_array2D_saving(i_file, Bphi, n_pts, n_psi, 'Bphi_mat'//char(0))

    ! 5. Close the file
    call HDF5_close(i_file)

#else
    write(*,*) 'Warning: HDF5 not enabled. File not saved.'
#endif
  end subroutine save_flux_data_h5


  subroutine save_avg_quantities_h5(filename, Psi, avg_vals)
        character(len=*), intent(in) :: filename
        real*8, intent(in) :: Psi(:)
        real*8, intent(in) :: avg_vals(:,:)

#ifdef USE_HDF5
        integer(HID_T) :: i_file
        integer :: ierr, n_psi

        ierr = 0
        n_psi = size(Psi)

        ! 1. Create the HDF5 file
        call HDF5_create(trim(filename), i_file, ierr)
        
        if (ierr /= 0) then
            write(*,*) 'Error: Could not create HDF5 file ', trim(filename)
            return
        end if

        ! 2. Save the radial coordinate
        call HDF5_array1D_saving(i_file, Psi, n_psi, 'Psi'//char(0))

        ! 3. Save each quantity in its own field
        ! column 1: Temperature
        call HDF5_array1D_saving(i_file, avg_vals(:,1), n_psi, 'T'//char(0))
        
        ! column 2: Electron Density
        call HDF5_array1D_saving(i_file, avg_vals(:,2), n_psi, 'ne'//char(0))
        
        ! column 3: Parallel Electric Field
        call HDF5_array1D_saving(i_file, avg_vals(:,3), n_psi, 'E_parallel'//char(0))

        ! 4. Close the file
        call HDF5_close(i_file)
#else
        write(*,*) 'Warning: HDF5 not enabled. Average data not saved.'
#endif
    end subroutine save_avg_quantities_h5

end module mod_save_flux_hdf5