! Filter one full-domain, cell-centred PadeOps binary budget field offline.
program dealias_field
   use mpi
   use iso_fortran_env, only: int64, error_unit
   use ieee_arithmetic, only: ieee_is_finite
   use kind_parameters, only: rkind
   use spectralMod, only: spectral
   use decomp_2d, only: decomp_info, decomp_2d_init, get_decomp_info, decomp_2d_finalize
   use decomp_2d_io, only: decomp_2d_read_one, decomp_2d_write_one
   implicit none

   type(spectral) :: spect
   type(decomp_info) :: grid
   real(rkind), allocatable :: field(:,:,:)
   complex(rkind), allocatable :: field_hat(:,:,:)
   character(len=4096) :: input_file='', output_file='', input_namelist='', iomsg=''
   integer :: nx=0, ny=0, nz=0, prow=0, pcol=0, dealiasType=0
   real(rkind) :: Lx=0, Ly=0, Lz=0, dealiasFact=2._rkind/3._rkind
   logical :: PeriodicInZ=.false., exists, periodicbcs(3)
   integer :: ierr, rank, nranks, unit, ios, arglen, nonfinite, any_nonfinite
   namelist /INPUT/ input_file, output_file, nx, ny, nz, Lx, Ly, Lz, &
                    prow, pcol, PeriodicInZ, dealiasFact, dealiasType

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
   if (command_argument_count() /= 1) call fail('Usage: dealias_field input_dealias.dat')
   call get_command_argument(1, input_namelist, length=arglen, status=ios)
   if (ios /= 0 .or. arglen > len(input_namelist)) call fail('Invalid input namelist path.')
   open(newunit=unit, file=trim(input_namelist), status='old', action='read', iostat=ios, iomsg=iomsg)
   if (ios /= 0) call fail('Opening namelist: '//trim(iomsg))
   read(unit, nml=INPUT, iostat=ios, iomsg=iomsg)
   if (ios /= 0) call fail('Reading INPUT namelist: '//trim(iomsg))
   close(unit)

   if (len_trim(input_file)==0 .or. len_trim(output_file)==0) call fail('Specify input_file and output_file.')
   if (trim(input_file)==trim(output_file)) call fail('Input and output must be different files.')
   if (min(nx,ny,nz)<2) call fail('nx, ny and nz must each be at least 2.')
   if (.not. all(ieee_is_finite([Lx,Ly,Lz,dealiasFact]))) call fail('Domain/filter values must be finite.')
   if (min(Lx,Ly,Lz)<=0) call fail('Lx, Ly and Lz must be positive.')
   if (dealiasFact<=0 .or. dealiasFact>1) call fail('dealiasFact must be in (0,1].')
   if (dealiasType/=0 .and. dealiasType/=1) call fail('dealiasType must be 0 (cutoff) or 1 (smooth).')
   if (PeriodicInZ .and. mod(nz,2)/=0) call fail('PeriodicInZ requires even nz.')
   if (prow<0 .or. pcol<0) call fail('prow and pcol must be nonnegative.')
   if ((prow==0) .neqv. (pcol==0)) call fail('Set both prow and pcol, or use 0 for both (automatic).')
   if (prow>0) then
      if (int(prow,int64)*pcol/=nranks) call fail('prow*pcol must equal the MPI rank count.')
   end if

   if (rank==0) then
      inquire(file=trim(input_file), exist=exists, iostat=ios)
      if (ios/=0 .or. .not. exists) call fail('Cannot access input_file: '//trim(input_file))
      ! STATUS=NEW also protects existing files reached through relative paths or links.
      open(newunit=unit, file=trim(output_file), status='new', access='stream', &
           form='unformatted', action='write', iostat=ios, iomsg=iomsg)
      if (ios/=0) call fail('Output must be a new file in an existing directory: '//trim(iomsg))
      close(unit)
      write(*,'(A)') 'Reading: '//trim(input_file)
      if (.not. PeriodicInZ .and. dealiasType==0) &
         write(*,'(A)') 'PadeOps horizontal cutoff is fixed at strict 2/3 (dealiasFact is ignored).'
      if (PeriodicInZ .and. dealiasType==0) &
         write(*,'(A)') 'PadeOps periodic cutoff uses dealiasFact*pi/dx in all three directions.'
      if (dealiasType==1) write(*,'(A)') 'PadeOps smooth filter ignores dealiasFact.'
   end if
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   periodicbcs = [.true., .true., PeriodicInZ]
   call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
   call get_decomp_info(grid)
   call spect%init('x', nx, ny, nz, Lx/real(nx,rkind), Ly/real(ny,rkind), Lz/real(nz,rkind), &
                   'FOUR', '2/3rd', dimTransform=2, fixOddball=.false., exhaustiveFFT=.false., &
                   init_periodicInZ=PeriodicInZ, dealiasF=dealiasFact, dealiasType=dealiasType)

   allocate(field(grid%xsz(1),grid%xsz(2),grid%xsz(3)))
   allocate(field_hat(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   call decomp_2d_read_one(1, field, trim(input_file), grid)
   nonfinite = 0
   if (.not. all(ieee_is_finite(field))) nonfinite = 1
   call MPI_Allreduce(nonfinite, any_nonfinite, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
   if (any_nonfinite/=0) call fail('Input contains NaN or infinity.')
   call spect%fft(field, field_hat)
   call spect%dealias(field_hat)
   call spect%ifft(field_hat, field)
   nonfinite = 0
   if (.not. all(ieee_is_finite(field))) nonfinite = 1
   call MPI_Allreduce(nonfinite, any_nonfinite, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
   if (any_nonfinite/=0) call fail('Filtering produced NaN or infinity.')
   call decomp_2d_write_one(1, field, trim(output_file), grid)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (rank==0) write(*,'(A)') 'Wrote filtered field: '//trim(output_file)
   deallocate(field, field_hat)
   call spect%destroy()
   call decomp_2d_finalize()
   call MPI_Finalize(ierr)

contains

   subroutine fail(reason)
      character(len=*), intent(in) :: reason
      integer :: abort_ierr
      write(error_unit,'(A,I0,2A)') 'dealias_field rank ', rank, ': ', trim(reason)
      call MPI_Abort(MPI_COMM_WORLD, 1, abort_ierr)
      stop 1
   end subroutine fail

end program dealias_field
