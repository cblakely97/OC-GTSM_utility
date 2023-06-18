!----------------------------------------------------------------------
!-----------------------------------------------------------------------
      PROGRAM OGCM_DL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      use datetime_module
      use MPI
      use netcdf
      use gsw_mod_toolbox, only: gsw_SA_from_SP, gsw_rho, gsw_CT_from_t&
          , gsw_Nsquared, gsw_mlp, gsw_geo_strf_dyn_height
      implicit none
!     This program downloads OGCM data (e.g. GOFS 3.1 HYCOM), and computes  
!     depth-averaged baroclinic pressure gradients, buoyancy frequencies, 
!     sea surface density, and mixed-layer depth ratio on the structured 
!     grid, outputting it to a compressed NetCDF file
!
!      -and (if required)- 
!
!     Computes depth-averaged momentum dispersion, momentum dispersion
!     as a quadratic bottom friction coefficient, and depth-averaged 
!     kinetic energy from OGCM 3D velocities

!     Use below two lines for output in single precision (saves space)
!      integer  :: nf90_type = nf90_float  ! Type to store matrix arrays in NetCDF
!      integer,parameter  :: sz = 4 !
!     Use below two lines for output in double precision
      integer  :: nf90_type = nf90_double ! Type to store matrix arrays in NetCDF
      integer,parameter  :: sz = 8 !

      character(len=80) :: BC3D_Name ! Name of the BC3D netCDF file 
      integer  :: BC3Ddir     ! orientaton of longitude 
      integer  :: NZ, NX, NY  ! dimensions of 3D BC model required
      integer  :: BC3D_DT = 3 ! Default delta t for the OGCM (hours)
      integer  :: dfl = 2     ! Compression level for NetCDF4
      integer  :: TMULT       ! BC3D_DT multiplier (to skip times)
      integer  :: is          ! Start of the x-dimension when writing out array
      integer  :: tsind = 1   ! Start time index in netCDF file
      integer  :: OutType = 1 ! Output type; = 1 baroclinicity, = 2 velocity
      character(len=16) :: TS, TE ! Start and end times
      character(len=3) :: BCServer = 'ftp' ! Server type
      real*8,parameter :: BPGlim = 0.1d0 !upper/lower bound for BPGs
      real*8,parameter :: LatUL = 89d0 !upper limit for Lat in haversine 
      real*8,allocatable,dimension(:) :: BC3D_Lon, BC3D_Lat, BC3D_Z
      real(sz),allocatable,dimension(:,:,:) :: BC3D_SP, BC3D_T,        &
                                               BC3D_BCP
      real(sz),allocatable,dimension(:,:) :: BC2D_NM, BC2D_NB, BC2D_BX,&
                                BC2D_SigTS, BC2D_BY, BC2D_MLD, BC2D_CD,&
                                BC2D_DispX, BC2D_DispY, BC2D_KE
      real*8,allocatable,dimension(:,:) :: DX, DY, DUU, DVV, DUV, DU,  &
                                           DV, B
!     Variables for temporal integration
      type(datetime) :: TSDT, TEDT, CurDT
      integer :: ierr, myProc, mnProc, nextProc, prevProc
      real*8,parameter :: RhoWat0 = 1d3, PI = 3.141592653589793D0,     &
                       deg2rad = PI/180d0, Rearth = 6378206.4d0,       &
                       G = 9.80665d0, rad2deg = 180d0/PI
      real(sz),parameter :: SigT0 = -999d0, DFV = 1d0 !DFV = 1d4 
      ! Netcdf output info
      character(len=40) :: BC2D_Name ! Name of Output File
      integer :: ncid, ncformat, timenc_id, timenc_dim_id, NX_dim_id,  &
              CD_id, NYYY_dim_id, NY_dim_id, NYY_dim_id, MLD_id, NM_id,&
              SigTS_id, BPGX_id, BPGY_id, NB_id, lon_id, lat_id, KE_id,&
              latc_id, lats_id, lonc_id, strlen_dim_id, DispX_id, &
              DispY_id, depth_id, ele_id, node_dim_id, nele_dim_id, &
              nvertex_dim_id, BCSL_id
!     CPB: Variables for outputting on ADCIRC grid. I did my best to
!     make sure that variable names match what is in ADCIRC source code.
      CHARACTER(len=40) :: fort14 ! name of fort.14 grid file
      INTEGER :: NP, NE
      REAL*8,ALLOCATABLE,DIMENSION(:)   :: SLAM, SFEA, DP
      REAL*8,PARAMETER :: SFEA0=0d0, SLAM0 = 0d0
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: NM, NeiTabEle
!     CPB: indices for spatial interpolation (use bilinear interpolation
!     routine from ADCIRC)
      INTEGER,ALLOCATABLE :: INDXY(:,:), INDZ(:)
      REAL*8,ALLOCATABLE  :: WEIGHTS(:,:)
!     Areas of elements and total area connected to each node
      REAL*8,ALLOCATABLE :: Areas(:), TotalArea(:)
!     edgelengths for use in calculating derivatives. This is made in
!     the Calc_Areas subroutine then used in Calc_Derivatives then it is
!     no longer needed and is deallocated
      REAL*8,ALLOCATABLE :: FDXE(:,:), FDYE(:,:), SFMXEle(:),&
                            SFMYEle(:), xve(:,:), yve(:,:)
!     Elemental derivatives of basis functions
      REAL*8,ALLOCATABLE,DIMENSION(:) :: Dphi1Dx, Dphi2Dx, Dphi3Dx, &
                                         Dphi1Dy, Dphi2Dy, Dphi3Dy
!     Intermediate values needed for calculating depth-averaged terms
!     namely baroclinic pressures, depth-resolved HYCOM velocities on 
!     ADCIRC grid and depth-averaged HYCOM velocities on adcirc grid
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: BCP_ADC, ADC3D_U, ADC3D_V
      REAL*8,ALLOCATABLE,DIMENSION(:) :: ADC2D_U, ADC2D_V
!     Depth integrated stuff for mom disp (Duu, Duv, Dvv)
      REAL*8,ALLOCATABLE,DIMENSION(:) :: DUU_ADC, DVV_ADC, DUV_ADC
!     Baroclinic pressure gradients for output
      REAL*8,ALLOCATABLE,DIMENSION(:) :: BPG_ADCx, BPG_ADCy
!     Momentum dispersion on adcirc grid for output
      REAL*8,ALLOCATABLE,DIMENSION(:) :: DispX_ADC, DispY_Adc
!     Other output variables
      REAL*8,ALLOCATABLE,DIMENSION(:) :: NB_ADC, NM_ADC, SigTS_ADC,&
                                         MLD_ADC, BCSL_ADC
      REAL*8 :: start,finish
!     Variables for NDCHL downloaded HYCOM files
!     derived type to hold info read in from the ogcm file
      TYPE :: OGCM_DATA
         INTEGER :: NSNAPS
         TYPE(DATETIME),ALLOCATABLE :: TSNAP(:)
         CHARACTER(LEN=40),ALLOCATABLE :: TSFILE(:)
         CHARACTER(LEN=40),ALLOCATABLE :: UVFILE(:)
         CHARACTER(LEN=40),ALLOCATABLE :: SSHFILE(:)
      END TYPE
      TYPE(OGCM_DATA) :: FN
      CHARACTER(LEN=16) :: OGCMFILE ! name of text file with OGCM info

!-----------------------------------------------------------------------

!.......Initialize MPI
      call MPI_Init(ierr)
      call MPI_Comm_Size (MPI_Comm_World,mnProc,ierr)   ! Get number of procs
      call MPI_Comm_Rank (MPI_Comm_World,myProc,ierr)   ! Get MPI rank
!      if (myProc.eq.0) then
!         write(6,*) 'mnProc = ',mnProc
!      endif
!      write(6,*) 'myProc = ',myProc
!..... 
!     Get next and previous processors to inform when netCDF
!     file is free to write out into
      nextProc = myProc + 1
      if (nextProc.gt.mnProc-1) nextProc = 0    
      prevProc = myProc - 1
      if (prevProc.lt.0) prevProc = mnProc - 1    
 
!.......Read the Control File
      call Read_Input_File() 
!.....Read OGCM txt file with filenames (this is hardcoded for now)
      OGCMFILE='ogcm_data.txt'
      CALL READOGCMDATA(OGCMFILE)
!.......Start the download and process loop
      call OGCM_Run()
  
!.......Finish up the program 
      call MPI_Finalize(ierr)
      stop
     
!-----------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ I N P U T _ F I L E
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Read_Input_File()
      implicit none
      logical :: Fexist
!
      ! Read from standard input
      if (myProc.eq.0) then
         read(5,*) ! Read blank
         read(5,'(A16)') TS;           
         read(5,'(A16)') TE;           
         TSDT = strptime(trim(TS),"%Y-%m-%d %H:%M")
         TEDT = strptime(trim(TE),"%Y-%m-%d %H:%M")
         write(6,*) TSDT%isoformat(' ')            
         write(6,*) TEDT%isoformat(' ')            
         ! Check this a valid time
         if (.not.TSDT%isvalid().or..not.TEDT%isvalid()) then
            write(6,*) 'ERROR: Invalid datetime TS or TE: '   &
                     //'Must be formatted as yyyy-MM-dd HH:mm'
            call MPI_Abort(MPI_Comm_World,0,ierr)
         endif
         read(5,*) TMULT; print *,'TMULT = ', TMULT          
         read(5,'(A3)') BCServer; print *, 'BCServer = ', BCServer
         read(5,*) OutType
         if (OutType == 0) then
            write(6,*) 'INFO: Compute & output buoyancy frequencies'
         elseif (OutType == 1) then
            write(6,*) 'INFO: Compute & output all T and S variables'
         elseif (OutType == 2) then
            write(6,*) 'INFO: Compute & output T, S, U, V'
         ELSEIF ( OutType.GE.3 ) THEN
            IF (OutType.EQ.3) THEN
               ! calculate ALL variables on ADCIRC grid
               WRITE(6,*) 'INFO: Compute & output T and S' &
                        // ' variables on ADCIRC grid'
            ELSEIF (OutType.EQ.4) THEN
               ! calculate ALL variables on ADCIRC grid
               WRITE(6,*) 'INFO: Compute & output T, S, U, and V' &
                        // ' variables on ADCIRC grid'
            ELSEIF (OutType.EQ.5) THEN
               ! calculate ALL variables on ADCIRC grid
               WRITE(6,*) 'INFO: Compute & output baroclinic sea' &
                       //' level on ADCIRC grid'
            ELSE
               write(6,*) 'ERROR: Invalid OutType = ',OutType
               call MPI_Finalize(ierr); stop
            ENDIF
            READ(5,'(A40)') fort14; print *, 'fort.14 = ',fort14
            ! read in fort.14
            CALL Read_f14()
            ! calculate both Areas and TotalArea
            CALL Calc_Areas()
            ! calculate dphi/dx, dphi/dy
            CALL Calc_Derivatives()
         else
            write(6,*) 'ERROR: Invalid OutType = ',OutType
            call MPI_Finalize(ierr); stop
         endif
         read(5,'(A40)') BC2D_Name; print *, 'BC2D_Name = ', BC2D_Name
         ncformat = nf90_hdf5
      endif
      ! Broadcast information to all processors
      call MPI_Bcast(TS,16,MPI_Character,0,MPI_Comm_World,ierr)
      call MPI_Bcast(TE,16,MPI_Character,0,MPI_Comm_World,ierr)
      TSDT = strptime(TS,"%Y-%m-%d %H:%M")
      TEDT = strptime(TE,"%Y-%m-%d %H:%M")
      call MPI_Bcast(TMULT,1,MPI_Integer,0,MPI_Comm_World,ierr)
      call MPI_Bcast(BCServer,3,MPI_Character,0,MPI_Comm_World,ierr)
      call MPI_Bcast(OutType,1,MPI_Integer,0,MPI_Comm_World,ierr)
      call MPI_Bcast(BC2D_Name,40,MPI_Character,0,MPI_Comm_World,ierr)
      IF (OutType.GE.3) THEN
         ! broadcast NP and NE
         CALL MPI_Bcast(NP,1,MPI_Integer,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(NE,1,MPI_Integer,0,MPI_Comm_World,ierr)
         ! allocate arrays on processors not equal to 0
         IF (myproc.gt.0) THEN
            ! allocate stuff
            ALLOCATE(SLAM(NP), SFEA(NP), DP(NP), NM(NE,3))
            ALLOCATE( Areas(NE), TotalArea(NP) )
            ALLOCATE( Dphi1Dx(NE), Dphi2Dx(NE), Dphi3Dx(NE), &
                      Dphi1Dy(NE), Dphi2Dy(NE), Dphi3Dy(NE) )
         ENDIF
         ! slam, sfea, DP
         CALL MPI_Bcast(SLAM,NP,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(SFEA,NP,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(DP,NP,MPI_REAL,0,MPI_Comm_World,ierr)
         ! NM
         CALL MPI_Bcast(NM,size(NM),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         ! Areas, TotalArea
         CALL MPI_Bcast(Areas,NE,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(TotalArea,NP,MPI_REAL,0,MPI_Comm_World,ierr)
         ! dphi_i
         CALL MPI_Bcast(Dphi1Dx,NE,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(Dphi2Dx,NE,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(Dphi3Dx,NE,MPI_REAL,0,MPI_Comm_World,ierr) 
         CALL MPI_Bcast(Dphi1Dy,NE,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(Dphi2Dy,NE,MPI_REAL,0,MPI_Comm_World,ierr)
         CALL MPI_Bcast(Dphi3Dy,NE,MPI_REAL,0,MPI_Comm_World,ierr) 
      ENDIF
!
      end subroutine Read_Input_File
!
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ I N P U T _ F I L E
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine OGCM_Run()
      
      implicit none
      type(timedelta) :: EndCheck
      character(len=3):: hhh 
      integer :: IT, ii
      logical :: OKflag

      ! Set the output filename
      write(hhh,'(I3)') myProc
      BC3D_Name = 'HYCOM'//trim(adjustl(hhh))//'.nc' 
      CurDT = TSDT + timedelta(hours=myProc*BC3D_DT*TMULT)
      write(6,*) 'MyProc = ',myProc,'CurDT = ',CurDT%isoformat(' ')
      EndCheck = TEDT - CurDT; IT = 0
      do while (EndCheck%total_seconds() >= 0)
         !call cpu_time(start)
         IF (OutType.LT.3) THEN
            ! output on native hycom grid
            do ii = 1,max(OutType,1)
               ! Download new OGCM NetCDF file
               call BC3DDOWNLOAD(CurDT,ii,OKflag)            
               if (.not.OKflag) exit
               ! Read the OGCM NetCDF file
               call Read_BC3D_NetCDF(ii)
               ! Calculate the new BC2D terms.
               call Calc_BC2D_Terms(ii)
            end do
         ELSEIF (OutType.EQ.5) THEN
            CALL GETOGCMFILE(1,CurDT,BC3D_Name)
            OKflag = .true.
            WRITE(6,*) 'BC3D_NAME = ',trim(BC3D_NAME)
            BC3D_NAME = trim(BC3D_NAME)
            ! Read the OGCM NetCDF file
            call Read_BC3D_NetCDF(3)
            ! Calculate the steric adjustment.
            call Calc_Steric_Adjustment()
            WRITE(6,*) "Calced BCSL for this time step"
         ELSEIF (OutType.GE.3) THEN
            DO ii = 1,MAX(OutType-2,1) 
               ! output on adcirc grid
               ! Download new OGCM NetCDF file
               !call BC3DDOWNLOAD(CurDT,ii,OKflag)            
               !IF (.NOT.OKflag) then
               !   WRITE(6,*) 'Issue with downloading '
               !   IF (IT.EQ.0) THEN
               !      WRITE(6,*) 'First timestep, can''t work'
               !   ELSE
               !      WRITE(6,*) 'Will use previous timestep data'
               !      Okflag = .TRUE.
               !   ENDIF
               !   exit
               !ENDIF
               ! get the filename from our txt file from earlier (this
               ! is also hardcoded for now)
               CALL GETOGCMFILE(ii,CurDT,BC3D_Name)
               OKflag = .true.
               WRITE(6,*) 'BC3D_NAME = ',trim(BC3D_NAME)
               BC3D_NAME = trim(BC3D_NAME)
               ! Read the OGCM NetCDF file
               call Read_BC3D_NetCDF(ii + 2)
               ! Calculate the new BC2D terms.
               call Calc_BC2D_Terms(ii + 2)
            ENDDO
         ENDIF
         ! Put new BC2D terms in netCDF output file
         call UpdateNetCDF(IT,OKflag)
         !call cpu_time(finish)
         !write(6,*) 'elapsed time = ',finish-start,'seconds'
         ! Update the time 
         CurDT = CurDT + timedelta(hours=mnProc*BC3D_DT*TMULT)
         EndCheck = TEDT - CurDT
         IT = IT + 1
      enddo
!
      end subroutine OGCM_Run
!
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ B C 3 D _ N E T C D F
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Read_BC3D_NetCDF(flag)
      implicit none
      integer,intent(in) :: flag
      logical :: errL
      integer :: NC_ID
      
      ! Open NETCDF 
      call Check_err(NF90_OPEN(BC3D_Name,NF90_NOWRITE,NC_ID))
    
      ! Get the dimensions and time and spatial arrays of the data
      ! or check the newly downloaded data
      call Get_LonLatDepthTime(NC_ID)
     
      ! Read all the necessary temporal varying data
      select case(flag)
        case(1) ! Baroclinicity
          write(6,*) myProc,'Processing T & S'
          ! Practical Salinity
          BC3D_SP = read_nc_var(NC_ID,'salinity  ',1)
          ! Temperature 
          BC3D_T  = read_nc_var(NC_ID,'water_temp',1)
        case(2) ! Velocity
          write(6,*) myProc,'Processing U & V'
          ! East-West velocity
          BC3D_SP = read_nc_var(NC_ID,'water_u   ',1)
          ! North-South Velocity
          BC3D_T  = read_nc_var(NC_ID,'water_v   ',1)
        CASE(3)
          write(6,*) myProc,'Processing T & S (ADCIRC grid)'
          ! Practical Salinity
          BC3D_SP = read_nc_var(NC_ID,'salinity  ',1)
          ! Temperature 
          BC3D_T  = read_nc_var(NC_ID,'water_temp',1)
        case(4) ! Velocity
          write(6,*) myProc,'Processing U & V (ADCIRC grid)'
          ! East-West velocity
          BC3D_SP = read_nc_var(NC_ID,'water_u   ',1)
          ! North-South Velocity
          BC3D_T  = read_nc_var(NC_ID,'water_v   ',1)
      end select

      ! Close NETCDF
      call Check_err(NF90_CLOSE(NC_ID))
!
!-----------------------------------------------------------------------
      end subroutine Read_BC3D_NetCDF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Get_LonLatDepthTime(NC_ID)
      implicit none
      integer,intent(in) :: NC_ID
      integer  :: Temp_ID, i, j, ii 
      real*8   :: BC3D_Lon_s, LatVal 
      real*8,allocatable :: templon(:)
      LOGICAL,SAVE :: FirstCall = .TRUE.
      INTEGER :: NY_temp, NX_temp ! to see if we swapped grids
        
      ! if this is not the first call, check to see if we have crossed
      ! over the time period from the GLBv0.08 grid to the GLBy0.08
      ! grid. The way I will check this is if NY has changed (since only
      ! resolution in the latitudinal direction changed). I am also gonna
      ! code it such that if it is trying to output on the HYCOM grid
      ! and we switch grids it aborts the program and writes an error
      ! message. If we are using the ADCIRC grid then it will recompute
      ! the interpolants and continue outputting.
      IF (.NOT.FirstCall) THEN
         !WRITE(6,*) 'checking grid'
         call Check_err(NF90_INQ_DIMID(NC_ID,'lat',Temp_ID))
         !write(6,*) 'inquiring dim'
         Call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,&
                        len=NY_temp))
         IF (NY_temp.NE.NY.AND.OutType.LT.3) THEN
            WRITE(6,*) 'ERROR: This time period contains the switch'&
                        //'from GLBv0.08 to GLBy0.08 grids.'
            WRITE(6,*) '       Currently there is no support for'&
                          //' OutType 1 or 2 for this swap.'
            CALL MPI_Finalize(ierr);stop
         ELSEIF (NY_temp.NE.NY) THEN
            ! deallocate all values we get from the input data
            write(6,*) 'INFO: we have switched grids, recalculating'&
                       //' interpolation weights.'
            DEALLOCATE( BC3D_Lat, BC3D_Lon, BC3D_Z, BC3D_SP, BC3D_T )
            FirstCall = .True.
         ENDIF
         ! test if we are still in the same lon coordinates (-180/180
         ! or 0/360)
         write(6,*) 'inquiring lon name'
         call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
         write(6,*) 'inquiring lon dimension'
         !Call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,&
         !               len=NY_temp))
         write(6,*) 'NX_temp = ',ny_temp
         write(6,*) 'NX = ',NX
         ALLOCATE( templon(NX) )
         call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
         !write(6,*) 'reading in templon'
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,templon))
         write(6,*) 'templon: ',minval(templon)
         IF (ALLOCATED(BC3D_LON)) THEN
            write(6,*) 'bc3d_lon: ',minval(BC3D_Lon)
            IF (minval(templon).ne.minval(BC3D_Lon)) THEN
               write(6,*) 'INFO: we have switched grids, recalculating'&
                          //' interpolation weights.'
               ! deallocate all values we get from the input data
               DEALLOCATE( BC3D_Lat, BC3D_Lon, BC3D_Z, BC3D_SP, BC3D_T )
               FirstCall = .True.
            ENDIF
         ENDIF
         DEALLOCATE( templon )
      ENDIF
      ! if this is the first call or we swap grids define all the
      ! various arrays and sizes we need
      IF (FirstCall) THEN
         FirstCall = .FALSE.
         ! First call
         call Check_err(NF90_INQ_DIMID(NC_ID,'lat',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
         call Check_err(NF90_INQ_DIMID(NC_ID,'lon',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))
         call Check_err(NF90_INQ_DIMID(NC_ID,'depth',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NZ))
          
         ! Allocate the lat lon and z arrays first 
         allocate(BC3D_Lat(NY),BC3D_Lon(NX),BC3D_Z(NZ))
         allocate(BC3D_SP(NX,NY,NZ),BC3D_T(NX,NY,NZ))
         ! only need to allocate BC2D variables if OutType < 3
         IF (OutType.LT.3) THEN
            ALLOCATE( BC3D_BCP(NX,NY,NZ) )
            allocate(BC2D_NB(NX,NY),BC2D_NM(NX,NY),&
                     BC2D_SigTS(NX,NY),BC2D_MLD(NX,NY))
            allocate(BC2D_BX(NX,NY),BC2D_BY(NX,NY-1),DX(NX,NY),&
                     DY(NX,NY-1))
         ENDIF
         if (OutType.eq.2) then
            allocate(BC2D_CD(NX,NY-2),BC2D_DispX(NX,NY-2),&
                     BC2D_DispY(NX,NY-2))
            allocate(DUU(NX,NY), DVV(NX,NY), BC2D_KE(NX,NY),&
                     DUV(NX,NY), DU(NX,NY), DV(NX,NY), B(NX,NY))
         endif

         ! Read the latitude, longitude and z variables 
         call Check_err(NF90_INQ_VARID(NC_ID,'lat',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lat))
         call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lon))
         call Check_err(NF90_INQ_VARID(NC_ID,'depth',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Z, &
                        start=[1],count=[NZ]))
         IF (OutType.LT.3) THEN
         ! Pre-computing the distances on the sphere
            do j = 1,NY 
               do i = 1,NX
                  LatVal  = min(LatUL,BC3D_Lat(j))
                  ii = i + 1
                  if (i.eq.NX) ii = 1
                  DX(i,j) = haversine(BC3D_Lon(i),BC3D_Lon(ii),   &
                                      LatVal,LatVal) 
                  if (j < NY) then
                     DY(i,j) = haversine(BC3D_Lon(i),BC3D_Lon(i), &
                                         BC3D_Lat(j),BC3D_Lat(j+1)) 
                  endif 
               enddo
            enddo
         ENDIF
         ! If we are outputting to ADCIRC grid just put this all in
         ! 0-360 coordinates 
         IF (OutType.GE.3) THEN
            DO i = 1,NY
               !IF (BC3D_Lon(i).lt.0d0) THEN
               !   BC3D_Lon(i) = BC3D_Lon(i) + 360d0
               !ENDIF
            ENDDO
            ! get weights and indices for linear interpolation
            WRITE(6,*) 'INFO: calculating interpolants'
            CALL Get_Interpolation_Weights()
         ENDIF
      endif
      BC3D_Lon_s = BC3D_Lon(1)
      if (BC3D_Lon_s < 0d0) then
         ! Represents seam at 180/-180
         BC3Ddir = -1
         !write(6,*) 'Lon of ',trim(BC3D_Name),&
         !           ' is defined from -180 to 180'
         is = 1
      else
         ! Represents seam at 0/360
         BC3Ddir =  1
         is = maxloc(BC3D_Lon,dim=1,mask=BC3D_Lon < 180d0) + 1
         !write(6,*) 'Lon of ',trim(BC3D_Name),&
         !           ' is defined from 0 to 360, is = ',is
      endif
      
      
      end subroutine Get_LonLatDepthTime
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      function read_nc_var(NC_ID,varname,BCIT_IN) result (Var)
      implicit none

      integer,intent(in) :: NC_ID, BCIT_IN
      character(10),intent(in) :: varname
      integer :: Temp_ID, i, j, k 
      real(sz) :: FV, SF, OS
      real(sz),allocatable :: Var(:,:,:)

      call Check_err(NF90_INQ_VARID(NC_ID,trim(varname),Temp_ID))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'_FillValue',FV))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'scale_factor',SF))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'add_offset',OS))
      allocate(Var(NX,NY,NZ)) 
      call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var))
      ! Add on a little for buffer
      FV = FV + 1d-3
      do j = 1,NY
         do i = 1,NX
            do k = 1,NZ
               if (Var(i,j,k) > FV) Var(i,j,k) = Var(i,j,k)*SF+OS
            enddo
         enddo
      enddo
!
!-----------------------------------------------------------------------
      end function read_nc_var
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E   C H E C K  _  E R R
!-----------------------------------------------------------------------
!     jgf49.17.02 Checks the return value from netCDF calls; if there
!     was an error, it writes the error message to the screen
!-----------------------------------------------------------------------
      subroutine check_err(iret)
      implicit none
      REAL*8,ALLOCATABLE,DIMENSION(:) :: test
      integer, intent(in) :: iret

      if (iret .ne. nf90_noerr) then
         write(6,*) 'netcdf check_err error' 
         test(1) = 0d0
         call MPI_Abort(MPI_Comm_World,iret,ierr)
      endif
!-----------------------------------------------------------------------
      end subroutine check_err
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine BC3DDOWNLOAD(TimeIn,flag,OKflag)
      use ifport, only: systemqq, sleep, getlasterrorqq
      
      implicit none
     
      type(datetime),intent(in) :: TimeIn
      integer,intent(in) :: flag
      logical,intent(out) :: OKflag
      logical(4) :: resOK
      type(datetime)     :: TimeOff
      character(len=200) :: fileserver, vars, times, filemid, & 
                            fileend, expt, options
      character(len=280) :: FullPath, FullPathP, line
      character(len=3)   :: hhh 
      character(len=4)   :: yyyy 
      character(len=5)   :: command 
      integer            :: iter, stat, fid, resINT, &
                            mpistat(mpi_status_size)
      integer*8          :: Fsize, FSizeTrue,file_exists
       
      ! OKflag is false by default    
      OKflag = .false.
      
      ! Add on the time based on processor number 
      TimeOff  = TimeIn - timedelta(hours=12)
      yyyy     = TimeOff%strftime("%Y")

      ! GLBv0.08: GOFS 3.1 on lat lon structured grid
      ! CPB 1/21/2022: changed to updated GOFS 3.1 grid (GLBy0.08)
      ! Get servername
      if (BCServer == 'ftp') then
         ! FTP
         ! As of February 18, 2020 the GLBv0.08 grid was discontinued
         ! and a new grid was implemented. Therefore we need to check to
         ! see if that is the case here. 
         ! NOTE: I have only implemented this update for FTP downloads!
         ! 
         ! To simplify the logic just assume we can use the GLBv0.08
         ! grid and then check if we need to switch to GLBy0.08 below
         fileserver = '"ftp://ftp.hycom.org/datasets/GLBv0.08/'
         IF ( TimeOff%getYear().GE.2021 ) THEN
            ! If we are in 2021 or later we must use GLBy0.08
            fileserver = '"ftp://ftp.hycom.org/datasets/GLBy0.08/'
         ELSEIF ( TimeOff%getYear().EQ.2020 ) THEN
            ! If we are in 2020 we need to check if we are after Feb 19.
            ! Start by checking if we are in March or later
            IF ( TimeOff%getMonth().GE.3 ) THEN
               ! If we are in March or later just use GLBy0.08
               fileserver = '"ftp://ftp.hycom.org/datasets/GLBy0.08/'
            ELSEIF (TimeOff%getMonth().EQ.2) THEN
               ! If we are in February check if we are on or after the
               ! 19th
               IF ( TimeOff%getDay().GE.19 ) THEN
                  ! we need to use GLBy0.08 grid
                  fileserver = '"ftp://ftp.hycom.org/datasets/GLBy0.08/'
               ENDIF
            ENDIF
         ENDIF
         command    = 'curl '
         options    = ' -s --connect-timeout 30'
      elseif (BCServer == 'ncs') then
         ! NCSS
         fileserver = '"http://ncss.hycom.org/thredds/ncss/GLBv0.08/'
         command    = 'curl '
         options    = ' -s --connect-timeout 30'
      elseif (BCServer == 'opd') then
         ! OpenDap
         command    = 'ncks '
         fileserver = '"http://tds.hycom.org/thredds/dodsC/GLBv0.08/'
         options    = ' -4 -O '
      endif
      if (TimeOff%getYear().ge.2018) then
         ! GOFS 3.1 analysis 
         expt = 'expt_93.0'
         filemid = '/hycom_glbv_930_'
         ! CPB: use same logic above to swap to glby if we are after
         ! February 19, 2020
         IF ( TimeOff%getYear().GE.2021 ) THEN
            filemid = '/hycom_glby_930_'
         ELSEIF ( TimeOff%getYear().EQ.2020 ) THEN
            IF (TimeOff%getMonth().GE.3 ) THEN
               filemid = '/hycom_glby_930_'
            ELSEIF (TimeOff%getMonth().EQ.2 ) THEN
               IF ( TimeOff%getDay().GE.19 ) THEN
                  filemid = '/hycom_glby_930_'
               ENDIF
            ENDIF
         ENDIF
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2017) then
         ! GOFS 3.1 analysis
         if (TimeOff%getMonth().ge.10) then
            expt = 'expt_92.9'
            filemid = '/hycom_glbv_929_'
         elseif (TimeOff%getMonth().ge.6) then
            expt = 'expt_57.7'
            filemid = '/hycom_GLBv0.08_577_'
         elseif (TimeOff%getMonth().ge.2) then
            expt = 'expt_92.8'
            filemid = '/hycom_glbv_928_'
         else
            expt = 'expt_57.2'
            filemid = '/hycom_GLBv0.08_572_'
         endif
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2016.and.TimeOff%getMonth().ge.5) then
         ! GOFS 3.1 analysis
         filemid = '/hycom_GLBv0.08_572_'
         expt = 'expt_57.2'
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2016.or.(TimeOff%getYear().eq.2015.and.&
              TimeOff%getYear().eq.12.and.TimeOff%getDay().eq.31)) then
         ! GOFS 3.1 analysis
         filemid = '/hycom_GLBv0.08_563_'
         expt = 'expt_56.3'
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().ge.1994.or.&
             (TimeOff%getYear().eq.1993.and.TimeOff%getMonth().ge.10)) then
         ! GOFS 3.1 reanalysis, expt 53.X with the year
         expt = 'expt_53.X/data/'//yyyy
         if (TimeOff%getYear().lt.2000) then
            filemid = '/hycom_GLBv0.08_530_'
         elseif (TimeOff%getYear().eq.2000) then
            filemid = '/hycom_GLBv0.08_531_'
         elseif (TimeOff%getYear().lt.2004) then
            filemid = '/hycom_GLBv0.08_532_'
         elseif (TimeOff%getYear().lt.2006.and.TimeOff%getMonth().lt.7) then
            filemid = '/hycom_GLBv0.08_533_'
         elseif (TimeOff%getYear().lt.2008) then
            filemid = '/hycom_GLBv0.08_534_'
         elseif (TimeOff%getYear().lt.2010) then
            filemid = '/hycom_GLBv0.08_535_'
         elseif (TimeOff%getYear().lt.2012) then
            filemid = '/hycom_GLBv0.08_536_'
         elseif (TimeOff%getYear().lt.2014) then
            filemid = '/hycom_GLBv0.08_537_'
         elseif (TimeOff%getYear().eq.2014) then
            filemid = '/hycom_GLBv0.08_538_'
         elseif (TimeOff%getYear().lt.2015) then
            filemid = '/hycom_GLBv0.08_539_'
         endif
      else
         write(6,*) 'Error: no GOFS 3.1 data available before Oct 1993.'
      endif
      hhh = ''
      write(hhh,'(I3)') myProc
      if (BCServer == 'ftp') then
         vars    = ''
         times   = TimeOff%strftime("%Y%m%d")//'12_t0'&
                 //TimeOff%strftime("%H")
         if (expt(6:6) == '9') then
            if (flag == 1) then
               fileend = '_ts3z.nc"'
            elseif (flag == 2) then
               fileend = '_uv3z.nc"'
            endif
         else
            if (flag == 2) then
               OKflag = .true.; return
            endif 
            fileend = '.nc"'
         endif
      elseif (BCServer == 'ncs') then
         if (flag == 1) then
            vars = '?var=salinity&var=water_temp'
         elseif (flag == 2) then
            vars = '?var=water_u&var=water_v'
         endif
         times   = '&time='//TimeIn%strftime("%Y-%m-%dT%H")&
                 //'%3A00%3A00Z'
         filemid = ''
         fileend = '&vertStride=1&addLatLon=true&accept=netcdf4"'
      elseif (BCServer == 'opd') then
         if (flag == 1) then
            vars = '" -v salinity,water_temp'
         elseif (flag == 2) then
            vars = '" -v water_u,water_v'
         endif
         times   = ' -d time,'//trim(hhh)
         fileend = ''
         filemid = ''
      endif
      ! Get the final path
      Fullpath = trim(fileserver)//trim(expt)//trim(vars)       &
               //trim(filemid)//trim(times)//trim(fileend) 
      WRITE(6,*) Fullpath
      ! Let's get expected filesize
      if (BCServer == 'ftp') then
         FSizeTrue = 0; 
         !do while(FSizeTrue.eq.0) 
            ! Get filesize at remote location 
            resOK = systemqq(command//trim(Fullpath)//trim(options)//  &
                       ' -I -L -o filesize'//trim(adjustl(hhh))//'.txt')
            if (.not.resOK) then
               resINT = getlasterrorqq()
               write(6,*) 'Error: ',resINT,myProc
            endif
            open(newunit=fid,                                     &
                 file='filesize'//trim(adjustl(hhh))//'.txt',     &
                 status='old',action='read',IOSTAT=file_exists)
               do ! read until eof
                  read(fid,'(A280)',iostat=stat) line
                  if (stat.ne.0) exit
                  if (line(1:15).eq.'Content-Length:') then
                     read(line(17:280),*) FSizeTrue
                     exit
                  endif
               enddo
            close(fid,STATUS='DELETE')
            if (FSizeTrue.eq.0) then
               ! Non-existant file. 
               ! Try forecast
               !TimeOff = TimeOff - timedelta(hours=24)
               !times   = TimeOff%strftime("%Y%m%d")//'12_t0'       &
               !        //TimeOff%getHour + 24
               !Fullpath = trim(fileserver)//trim(expt)//trim(vars) &
               !        //trim(filemid)//trim(times)//trim(fileend)
               !if (FullPath.eq.FullPathP) then
               !   write(6,*) myProc,                               &
               !             'Previous and current paths are equal'
               !endif
               ! Try BC3D_DT hours ahead
               !TimeOff = TimeOff + timedelta(hours=BC3D_DT)
               !times   = TimeOff%strftime("%Y%m%d")//'12_t0'       &
               !        //TimeOff%strftime("%H")
               !Fullpath = trim(fileserver)//trim(expt)//trim(vars) &
               !        //trim(filemid)//trim(times)//trim(fileend)
               write(6,*) myProc,'no data available at: ',trim(FullPath)
               return ! return error
            endif 
         !enddo
      else
         ! Expect over a 500MB
         FSizeTrue = 5e8
      endif
      
      write(6,*) myProc, 'Downloading: ',trim(FullPath)
      iter = 0; resOK = .false. 
      do while (.not.resOK.and.iter < 25) 
         resOK = systemqq(command//trim(Fullpath)//trim(options)//&
                         ' -o'//trim(BC3D_Name))
         if (resOK) then
            ! Expect over a 500MB
            inquire(file=trim(BC3D_Name),size=Fsize)
            if (Fsize < FSizeTrue) then
               write(6,*) 'FSize of '//trim(BC3D_Name)//': ',     &
                           FSize, FSizeTrue
               resOK = .false.
            endif
         endif
         if (.not.resOK) then
            ! Let's wait 5 seconds before trying again
            write(6,*) 'Problem downloading GOFS 3.1 NetCDF. '&
                    // 'Try again after 5 s'
            call sleep(5) ; iter = iter + 1; 
         endif
      enddo
      if (.not.resOK) then
         write(6,*) 'Error in downloading GOFS 3.1 NetCDF.'
         call MPI_Abort(MPI_Comm_World,0,ierr) 
      else
         write(6,*) myProc,'download successfull'
         OKflag = .true.
      endif
!----------------------------------------------------------------------
      end subroutine BC3DDOWNLOAD
!-----------------------------------------------------------------------
!     S U B R O U T I N E  C A L C _ B C 2 D _ T E R M S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Calc_BC2D_Terms(flag)
      implicit none
     
      integer,intent(in) :: flag 
      real*8 :: rho, rhoavg, bx_avg, by_avg, bx_z, by_z, DZ
      real*8 :: uavg, vavg, DUUX, DVVY, DUVX, DUVY, Vel, CD, CDx, CDy 
      real*8 :: SA(NZ), CT(NZ), Lat(NZ), N2(NZ-1), zmid(NZ-1)
      integer  :: i, j, k, ii, kb, kbx, kby, ip, im
      real(sz) :: FV = -3d4, FVP = -3d4 + 1d-3, Vs = 1d-3 

      select case(flag)
      
      !! FOR SALINITY & TEMPERATURE 
      case(1)
      ! Loop over all cells to get the baroclinic pressure and the
      ! buoyancy frequencies
!!$omp parallel private(Lat,kb,DZ,rho,rhoavg,kbx,kby,bx_avg,by_avg, &
!!$omp bx_z,by_z,SA,CT,N2,zmid)
!!$omp do
      do j = 1,NY
         ! Need an array of latitude values for Nsquared call
         Lat = BC3D_Lat(j)
         do i = 1,NX
            ! Initialize the free surface density
            BC2D_SigTS(i,j) = SigT0
            kb = 0; 
            do k = 1,NZ
               if (k > 1) DZ = BC3D_Z(k) - BC3D_Z(k-1)
               ! For the baroclinic pressure gradient and surface
               ! density 
               if (BC3D_SP(i,j,k) > FVP .and. BC3D_T(i,j,k) > FVP) then
                  ! Change SP to SA
                  SA(k) = gsw_SA_from_SP(max(2d0,BC3D_SP(i,j,k)),&
                          BC3D_Z(k),BC3D_Lon(i),BC3D_Lat(j))
                  ! Change T to CT
                  CT(k) = gsw_CT_from_T(SA(k),dble(BC3D_T(i,j,k)),&
                          BC3D_Z(k))
                  ! Save previous rho
                  if (k > 1) rhoavg = rho
                  ! Calculate rho
                  rho = gsw_rho(SA(k),CT(k),BC3D_Z(k))
                  if (k > 1) then
                     ! The trapezoidal rule
                     rhoavg = 0.5d0*(rhoavg + rho)
                     BC3D_BCP(i,j,k) = BC3D_BCP(i,j,k-1)    &
                                     + DZ*(rhoavg - RhoWat0)
                  else
                     BC2D_SigTS(i,j) = rho - RhoWat0
                     BC3D_BCP(i,j,k) = 0.0d0
                  endif
                  kb = k
               else
                  ! Reached the bottom
                  BC3D_BCP(i,j,k) = FV
               endif
            enddo
            ! The buoyancy frequencies
            if (kb > 1) then
               ! Get the buoyancy frequency profile
               call gsw_Nsquared(SA(1:kb),CT(1:kb),BC3D_Z(1:kb),  &
                                 Lat(1:kb),N2(1:kb-1),zmid(1:kb-1))
               ! Get bottom (seabed) value
               BC2D_NB(i,j) = sqrt(max(0.0d0,N2(kb-1)))
               ! Integrate (here N2 is defined in the middle of depth
               ! layers (zmid), and I assume it to be the average in
               ! that depth range for simplicity) 
               BC2D_NM(i,j)  = 0.0d0  
               do k = 1,kb-1
                  DZ = BC3D_Z(k+1) - BC3D_Z(k)
                  BC2D_NM(i,j) = BC2D_NM(i,j) + DZ*sqrt(max(0.d0,N2(k)))
               enddo
               ! Divide by the depth
               BC2D_NM(i,j) = BC2D_NM(i,j)/BC3D_Z(kb)
               ! Get the mixed-layer depth
               !BC2D_MLD(i,j) = min(BC3D_Z(kb),&
               !                gsw_mlp(SA(1:kb),CT(1:kb),BC3D_Z(1:kb)))
               ! Get the mixed-layer depth ratio
               BC2D_MLD(i,j) = min(DFV,&
                     gsw_mlp(SA(1:kb),CT(1:kb),BC3D_Z(1:kb))/BC3D_Z(kb))
            else
               ! Let's just set to zero value 
               BC2D_NB(i,j)  = 0.0d0
               BC2D_NM(i,j)  = 0.0d0
               BC2D_MLD(i,j) = DFV
            endif
         enddo
      enddo 
!!$omp end do
      ! Now calculate the gradients of the BCP and integrate,
      ! and calculate gradients of DUU, DVU and DVV
      ! (we do central difference about the mid-point)
!!$omp do
      do j = 1,NY 
         do i = 1,NX
            kbx = 0; kby = 0;
            do k = 1,NZ
               if (k > 1) then
                  bx_avg = bx_z
                  by_avg = by_z 
                  ! If we hit bottom
                  if (BC3D_BCP(i,j,k) < FVP) exit
                  ii = i + 1
                  if (i == NX) ii = 1
                  if (BC3D_BCP(ii,j,k) > FVP) then
                     ! x-gradient at this level
                     bx_z = ( BC3D_BCP(ii,j,k) - BC3D_BCP(i,j,k) ) &
                          / DX(i,j)
                     ! trapezoidal rule
                     bx_avg = 0.5d0*(bx_avg + bx_z)
                     BC2D_BX(i,j) = BC2D_BX(i,j) +      &
                                   (BC3D_Z(k)-BC3D_Z(k-1))*bx_avg
                     kbx = k
                  endif
                  if (j < NY) then
                     if (BC3D_BCP(i,j+1,k) > FVP) then
                        ! y-gradient at this level
                        by_z = ( BC3D_BCP(i,j+1,k) -      &
                                 BC3D_BCP(i,j,k) ) / DY(i,j)
                        ! trapezoidal rule
                        by_avg = 0.5d0*(by_avg + by_z)
                        BC2D_BY(i,j) = BC2D_BY(i,j) +     & 
                                      (BC3D_Z(k)-BC3D_Z(k-1))*by_avg
                        kby = k
                     endif
                  endif
               else
                  ! Is always zero at top
                  bx_z = 0.0d0
                  by_z = 0.0d0
                  BC2D_BX(i,j) = 0.0d0
                  if (j < NY) BC2D_BY(i,j) = 0.0d0
               endif
            enddo
            ! Get the depth-averaged value
            if (kbx > 1) then
               BC2D_BX(i,j) = G/RhoWat0*min(BPGlim,max(-BPGlim, &
                              BC2D_BX(i,j)/BC3D_Z(kbx)))

            endif
            if (kby > 1) then
               BC2D_BY(i,j) = G/RhoWat0*min(BPGlim,max(-BPGlim, &
                              BC2D_BY(i,j)/BC3D_Z(kby)))
            endif
         enddo
      enddo
!!$omp end do
!!$omp end parallel

      !! FOR VELOCITIES 
      case(2)
!!$omp parallel private(kb,DZ,uavg,vavg,ip,im,DUUX,DUVY,DUVX,DVVY,Vel,CDx,CDy)
!!$omp do
      ! Loop over all cells to get the depth-integrated momentum dispersion
      do j = 1,NY
         do i = 1,NX
            ! Initialize the dispersion values
            DUU(i,j) = 0d0; DVV(i,j) = 0d0; DUV(i,j) = 0d0;
            DU(i,j) = 0d0; DV(i,j) = 0d0; B(i,j) = 0d0; 
            BC2D_KE(i,j) = 0d0; kb = 0; 
            do k = 2,NZ
               ! For depth-integrated dispersion 
               if (BC3D_SP(i,j,k) < FVP .or. BC3D_T(i,j,k) < FVP) exit
               DZ = BC3D_Z(k) - BC3D_Z(k-1)
               ! The trapezoidal rule
               uavg = 0.5d0*(BC3D_SP(i,j,k) + BC3D_SP(i,j,k-1))
               vavg = 0.5d0*(BC3D_T(i,j,k) + BC3D_T(i,j,k-1))
               DUU(i,j) = DUU(i,j) + DZ*uavg*uavg
               DVV(i,j) = DVV(i,j) + DZ*vavg*vavg
               DUV(i,j) = DUV(i,j) + DZ*uavg*vavg
               DU(i,j)  = DU(i,j) + DZ*uavg
               DV(i,j)  = DV(i,j) + DZ*vavg
               kb = k
            enddo
            ! Now calculate dispersion values in depth-averaged form
            if (kb > 1) then
               B(i,j)   = BC3D_Z(kb)
               ! Depth-averaged velocities
               DU(i,j)  = DU(i,j)/B(i,j)
               DV(i,j)  = DV(i,j)/B(i,j)
               ! Depth-averaged kinetic energy
               BC2D_KE(i,j)  = 0.5d0*( DUU(i,j) + DVV(i,j) ) / B(i,j)
               ! Depth-averaged dispersion
               DUU(i,j) = DUU(i,j)/B(i,j) - DU(i,j)*DU(i,j) 
               DVV(i,j) = DVV(i,j)/B(i,j) - DV(i,j)*DV(i,j)
               DUV(i,j) = DUV(i,j)/B(i,j) - DU(i,j)*DV(i,j)
            endif
         enddo
      enddo 
!!$omp end do

      ! Now calculate the gradients of DUU, DUV and DVV
      ! (we do central difference about original point)
!!$omp do
      do j = 2,NY-1 
         do i = 1,NX
            ip = i + 1
            im = i - 1
            if (ip > NX) ip = 1
            if (im < 1)  im = NX
                        
            if (abs(DU(i,j)) < Vs .or. abs(DV(i,j)) < Vs .or.      &
                B(ip,j) < BC3D_Z(2) .or. B(im,j) < BC3D_Z(2) .or.  & 
                B(i,j+1) < BC3D_Z(2) .or. B(i,j-1) < BC3D_Z(2)) then
               BC2D_CD(i,j-1) = 0d0
               BC2D_DispX(i,j-1) = 0d0
               BC2D_DispY(i,j-1) = 0d0
               cycle
            endif

            DUUX = 0.5d0*( DUU(ip,j) - DUU(im,j) ) / DX(i,j)
            DUVX = 0.5d0*( DUV(ip,j) - DUV(im,j) ) / DX(i,j)
            DVVY = 0.5d0*( DVV(i,j+1)  - DVV(i,j-1) ) / DY(i,j)
            DUVY = 0.5d0*( DUV(i,j+1)  - DUV(i,j-1) ) / DY(i,j)

            BC2D_DispX(i,j-1) = DUUX + DUVY
            BC2D_DispY(i,j-1) = DVVY + DUVX

            CD = BC2D_DispX(i,j-1)*DU(i,j) + BC2D_DispY(i,j-1)*DV(i,j)
            if (CD <= 0d0) then
               BC2D_CD(i,j-1) = 0d0
            else
               Vel  = sqrt(DU(i,j)*DU(i,j) + DV(i,j)*DV(i,j)) 
               BC2D_CD(i,j-1) = CD * B(i,j) / Vel**3
            endif 
            !CDx  = B(i,j) * ( DUUX + DUVY ) / ( DU(i,j) * Vel )
            !CDy  = B(i,j) * ( DVVY + DUVX ) / ( DV(i,j) * Vel )   
            !BC2D_CD(i,j-1) = sqrt(CDx*CDx + CDy*CDy)
            !BC2D_CD(i,j-1) = max(abs(CDx),abs(CDy))
         enddo
      enddo
!!$omp end do
!!$omp end parallel
      CASE(3)
      ! CPB: calculate T and S values on the ADCIRC mesh
         CALL Calc_BC2D_TS_ADCIRC()
      CASE(4)
         CALL Calc_BC2D_UV_ADCIRC()
      end select
!-----------------------------------------------------------------------
      end subroutine Calc_BC2D_Terms
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E  C A L C _ B C 2 D _ T S _ A D C I R C
!-----------------------------------------------------------------------
!     This subroutine takes the salinity and temperature read from GOFS
!     data and:
!        1) Converts from practical salinity to absolute salinity and
!           from temperature to conservative temperature on the GOFS
!           grid.
!        2) Loops through every ADCIRC node and for every depth in 
!           BC3D_Z less than the node depth. At each Z-level
!           interpolates the conservative temp and absolute salinity to
!           the grid point and calculates a baroclinic pressure.
!        3) Once the baroclinic pressures (BC3D_BCP) are calculated,
!           loop through once more to calculate the baroclinic pressure
!           gradients at each depth FOR EACH ELEMENT. We then integrate
!           over the depth and depth average using the shallowest depth
!           from the element (BC2D_BPGXe. 
!-----------------------------------------------------------------------
      SUBROUTINE Calc_BC2D_TS_ADCIRC()
      IMPLICIT NONE
      REAL*8 :: RHO, RHOAVG, DZ
      REAL*8 :: SA(NZ), CT(NZ), LAT(NZ), N2(NZ-1), ZMID(NZ-1)!,&
                !BCP_ADC(NZ,NP)
      REAL(SZ) :: FV = -3D4, FVP = -3D4 + 1D-3, VS = 1D-3 
      ! persistent logical for if we are on the first call
      LOGICAL,SAVE :: FirstCall = .TRUE.
      ! looping variables for do loops
      INTEGER :: IP, IE, iz, idzmax, truebottom, NM1, NM2, NM3
      ! for calculating baroclinic pressures
      REAL*8 :: xx, yy, sp1, sp2, sp3, sp4, t1, t2, t3, t4,&
                sa1, sa2, sa3, sa4, ct1, ct2, ct3, ct4
      ! for calculating baroclinic pressure gradients
      REAL*8 :: ddx1, ddx2, ddx3, ddy1, ddy2, ddy3, bpgx, bpgy,&
                bcp1, bcp2, bcp3, bx, by, facval
      ! allocate output arrays if on the first call
      IF (FirstCall) THEN
         ALLOCATE( BPG_ADCx(NP), BPG_ADCy(NP), SigTS_ADC(NP),&
                   NB_ADC(NP), NM_ADC(NP), MLD_ADC(NP) ) 
         FirstCall = .FALSE.
      ENDIF
      ALLOCATE( BCP_ADC(NZ,NP) )
      ! reset all values to 0
      BPG_ADCx  = 0d0
      BPG_ADCy  = 0d0
      BCP_ADC   = 0d0
      SigTS_ADC = 0d0
      NB_ADC    = 0d0
      NM_ADC    = 0d0
      MLD_ADC   = 0d0
      !
      ! get BCP_ADC, SigTS_ADC, NB_ADC, and NM_ADC
      !
      DO IP = 1,NP
         ! save lat and lon as yy,xx for readability
         xx = SLAM(IP)
         IF (minval(BC3D_Lon).GE.0d0.AND.xx.lt.0d0) xx = xx + 360d0

         yy = SFEA(IP)
         ! get max depth we are using
         idzmax = INDZ(IP)
         ! if we are only at the surface skip this whole loop
         IF (idzmax.lt.2) CYCLE
         ! reset the variable to track if we are at the bottom to 0
         truebottom = 0
         ! for my sanity reset the arrays for SA and CT to 0
         SA = 0d0
         CT = 0d0
         ! we also need a latitude array for the buoyancy frequency profile
         ! calculations
         Lat = yy
         DO iz = 1,idzmax
            ! for readability save the salinity and temps as individual values
            ! rather than using the indices
            ! salinity
            sp1 = BC3D_SP(INDXY(1,IP),INDXY(2,IP),iz)
            sp2 = BC3D_SP(INDXY(3,IP),INDXY(2,IP),iz)
            sp3 = BC3D_SP(INDXY(1,IP),INDXY(4,IP),iz)
            sp4 = BC3D_SP(INDXY(3,IP),INDXY(4,IP),iz)
            ! temperature
            t1  = BC3D_T(INDXY(1,IP),INDXY(2,IP),iz)
            t2  = BC3D_T(INDXY(3,IP),INDXY(2,IP),iz)
            t3  = BC3D_T(INDXY(1,IP),INDXY(4,IP),iz)
            t4  = BC3D_T(INDXY(3,IP),INDXY(4,IP),iz)
            ! if any of these are a fill value then we have hit the bottom and
            ! should use a fill value
            IF (sp1.lt.FVP .OR. sp2.lt.FVP .OR.&
                sp3.lt.FVP .OR. sp4.lt.FVP .OR.&
                t1.lt.FVP  .OR. t2.lt.FVP .OR.&
                t3.lt.FVP  .OR. t4.lt.FVP) THEN
               BCP_ADC(iz,IP) = FV
            ELSE
               ! convert to absolute salinity
               sa1 = gsw_SA_from_SP(max(2d0,sp1),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(1,IP)),BC3D_LAT(INDXY(2,IP)))
               !
               sa2 = gsw_SA_from_SP(max(2d0,sp2),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(3,IP)),BC3D_LAT(INDXY(2,IP)))
               !
               sa3 = gsw_SA_from_SP(max(2d0,sp3),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(1,IP)),BC3D_LAT(INDXY(4,IP)))
               !
               sa4 = gsw_SA_from_SP(max(2d0,sp4),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(3,IP)),BC3D_LAT(INDXY(4,IP)))
               ! convert to conservative temperature
               ct1 = gsw_CT_from_T(sa1,dble(t1),BC3D_Z(iz))
               ct2 = gsw_CT_from_T(sa2,dble(t2),BC3D_Z(iz))
               ct3 = gsw_CT_from_T(sa3,dble(t3),BC3D_Z(iz))
               ct4 = gsw_CT_from_T(sa4,dble(t4),BC3D_Z(iz))
               ! interpolate to our ADCIRC node
               SA(iz) = sa1*WEIGHTS(1,IP) + sa2*WEIGHTS(2,IP) +&
                        sa3*WEIGHTS(3,IP) + sa4*WEIGHTS(4,IP)
               CT(iz) = ct1*WEIGHTS(1,IP) + ct2*WEIGHTS(2,IP) +&
                        ct3*WEIGHTS(3,IP) + ct4*WEIGHTS(4,IP)
               ! if we are below the surface put rho from the last depth into
               ! rhoavg for trapezoidal rule later
               IF (iz.gt.1) rhoavg = rho
               ! get our new rho for this depth
               rho = gsw_rho(SA(iz),CT(iz),BC3D_Z(iz))
               ! if we are below the surface get baroclinic pressure,
               ! otherwise get surface values
               IF (iz.gt.1) THEN
                  rhoavg = 0.5d0*(rhoavg+rho)
                  dz = BC3D_Z(iz) - BC3D_Z(iz-1)
                  BCP_ADC(iz,IP) = BCP_ADC(iz-1,IP) +&
                                   dz*(rhoavg - RhoWat0)
               ELSE
                  SigTS_ADC(IP) = rho-RhoWat0
               ENDIF
               truebottom = iz
            ENDIF
         ENDDO
         ! if applicable get the buoyancy frequencies
         IF (truebottom.GT.1) THEN
            CALL gsw_Nsquared(SA(1:truebottom),CT(1:truebottom),&
                              BC3D_Z(1:truebottom),Lat(1:truebottom),&
                              N2(1:truebottom-1),zmid(1:truebottom-1))
            NB_ADC(IP) = SQRT(MAX(0d0,N2(truebottom-1)))
            ! integrate to get depth-averaged value. Since gsw_Nsquared
            ! returns the value at the midpoint I will assume that is
            ! the average in that z-level
            DO iz = 1,truebottom-1
               dz = BC3D_Z(iz+1)-BC3D_Z(iz)
               NM_ADC(IP) = NM_ADC(IP) + dz*sqrt(max(0d0,N2(iz)))
            ENDDO
            NM_ADC(IP) = NM_ADC(IP)/BC3D_Z(truebottom)
            ! Get Mixed-layer depth in case I decide I need it
            MLD_ADC(IP) = min(DFV,&
                             gsw_mlp(SA(1:truebottom),CT(1:truebottom),&
                               BC3D_Z(1:truebottom)))/BC3D_Z(truebottom)
         ENDIF
      ENDDO
      !
      ! Calculate baroclinic pressure gradients
      !
      DO IE = 1,NE
         NM1 = NM(IE,1)
         NM2 = NM(IE,2)
         NM3 = NM(IE,3)
         ddx1 = Dphi1Dx(IE)
         ddx2 = Dphi2Dx(IE)
         ddx3 = Dphi3Dx(IE)
         ddy1 = Dphi1Dy(IE)
         ddy2 = Dphi2Dy(IE)
         ddy3 = Dphi3Dy(IE)
         ! Only go down to the largest GOFS depth that is less than the ADCIRC
         ! depth
         idzmax = min(INDZ(NM1), INDZ(NM2), INDZ(NM3))-1
         ! reset bpgx and bpgy to 0
         bpgx = 0d0
         bpgy = 0d0
         bx = 0d0
         by = 0d0
         ! only calculate if we are deep enough
         IF (idzmax.gt.0) THEN
            DO iz = 1,idzmax
               bcp1 = BCP_ADC(iz+1,NM1)
               bcp2 = BCP_ADC(iz+1,NM2)
               bcp3 = BCP_ADC(iz+1,NM3)
               ! if we hit a fill value then stop this loop and fix the
               ! last value we added to the sum by subtracting off half
               ! of it
               IF (bcp1.lt.FVP.OR.bcp2.lt.FVP.OR.bcp3.lt.FVP) then
                  bpgx = bpgx - 0.5d0*bx*dz
                  bpgy = bpgy - 0.5d0*by*dz
                  idzmax = iz-1
                  EXIT
               ENDIF
               dz = BC3D_Z(iz+1)-BC3D_Z(iz)
               ! for ease of readability put bpg at this depth in variables
               ! called
               ! bx and by before adding them to the total sum
               bx = bcp1*ddx1 + bcp2*ddx2 + bcp3*ddx3
               by = bcp1*ddy1 + bcp2*ddy2 + bcp3*ddy3
               ! factor for trapezoidal rule
               facval = 1d0
               IF (iz.eq.idzmax) facval = 0.5d0 
               bpgx = bpgx + facval*bx*dz
               bpgy = bpgy + facval*by*dz
            ENDDO
            ! depth average
            IF (BC3D_Z(idzmax+1).GT.0) THEN
               bpgx = bpgx/BC3D_Z(idzmax+1)
               bpgy = bpgy/BC3D_Z(idzmax+1)
            ELSE
               bpgx = 0d0
               bpgy = 0d0
            ENDIF
         ENDIF
         ! area weight and sum
         BPG_ADCx(NM1) = BPG_ADCx(NM1) + Areas(IE)*bpgx
         BPG_ADCx(NM2) = BPG_ADCx(NM2) + Areas(IE)*bpgx
         BPG_ADCx(NM3) = BPG_ADCx(NM3) + Areas(IE)*bpgx
         BPG_ADCy(NM1) = BPG_ADCy(NM1) + Areas(IE)*bpgy
         BPG_ADCy(NM2) = BPG_ADCy(NM2) + Areas(IE)*bpgy
         BPG_ADCy(NM3) = BPG_ADCy(NM3) + Areas(IE)*bpgy
      ENDDO
      ! divide by totalarea
      DO IP = 1,NP
         BPG_ADCx(IP) = G/RhoWat0*BPG_ADCx(IP)
         BPG_ADCy(IP) = G/RhoWat0*BPG_ADCy(IP)
         IF (TotalArea(IP).GT.0d0) THEN
            BPG_ADCx(IP) = BPG_ADCx(IP)/TotalArea(IP)
            BPG_ADCy(IP) = BPG_ADCy(IP)/TotalArea(IP)
         ELSE
            BPG_ADCx(IP) = 0d0
            BPG_ADCy(IP) = 0d0
         ENDIF
      ENDDO 
      DEALLOCATE( BCP_ADC )
!-----------------------------------------------------------------------
      END SUBROUTINE Calc_BC2D_TS_ADCIRC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E  C A L C _ B C 2 D _ U V _ A D C I R C
!-----------------------------------------------------------------------
!     This subroutine takes the u-velocity and v-velocity read from GOFS
!     data and:
!        1) Interpolates the velocities to the ADCIRC grid in both 
!           depth-resolving and depth-averaged form
!        2) Calculates momentum dispersion on the ADCIRC grid in the
!           manner described in the ADCIRC theory and formulation report
!-----------------------------------------------------------------------
      SUBROUTINE Calc_BC2D_UV_ADCIRC()
      IMPLICIT NONE
      REAL*8 :: UU, UUAVG, VV, VVAVG, DZ
      REAL*8 :: Udiff, Udiffavg, Vdiff, Vdiffavg
      REAL(SZ) :: FV = -3D4, FVP = -3D4 + 1D-3, VS = 1D-3 
      ! persistent logical for if we are on the first call
      LOGICAL,SAVE :: FirstCall = .TRUE.
      ! looping variables for do loops
      INTEGER :: IP, IE, iz, idzmax, truebottom, NM1, NM2, NM3
      ! for interpolating velocities
      REAL*8 :: u1, u2, u3, u4, v1, v2, v3, v4
      ! for calculating momentum dispersion
      REAL*8 :: ddx1, ddx2, ddx3, ddy1, ddy2, ddy3, DUUx, DVVy,&
                DUVx, DUVy
      ! allocate output arrays if on the first call
      IF (FirstCall) THEN
         ALLOCATE( DispX_ADC(NP), DispY_ADC(NP) ) 
         FirstCall = .FALSE.
      ENDIF
      ! allocate what we need to the math in this subroutine
      ALLOCATE( ADC3D_U(NZ,NP), ADC3D_V(NZ,NP), &
                ADC2D_U(NP),    ADC2D_V(NP),    &
                DUU_ADC(NP),    DVV_ADC(NP),   &
                DUV_ADC(NP) )
      ! set errythign to 0
      DispX_ADC = 0d0
      DispY_ADC = 0d0
      ADC3D_U = 0d0
      ADC3D_V = 0d0
      ADC2D_U = 0d0
      ADC2D_V = 0d0
      DUU_ADC = 0d0
      DVV_ADC = 0d0
      DUV_ADC = 0d0
      !
      ! interpolate hycom velocities to adcirc grid and find
      ! depth-averaged hycom velocities on adcirc grid
      !
      DO IP = 1,NP
         ! find predetermined max depth index
         idzmax = INDZ(IP)
         ! if we are too shallow skip
         IF ( idzmax.LT.2 ) CYCLE
         ! in case there is a bathy mismatch
         truebottom = 0
         DO iz = 1,idzmax
            ! store individual values for readability
            ! u-velocity
            u1 = BC3D_SP(INDXY(1,IP),INDXY(2,IP),iz)
            u2 = BC3D_SP(INDXY(3,IP),INDXY(2,IP),iz)
            u3 = BC3D_SP(INDXY(1,IP),INDXY(4,IP),iz)
            u4 = BC3D_SP(INDXY(3,IP),INDXY(4,IP),iz)
            ! v-velocity
            v1  = BC3D_T(INDXY(1,IP),INDXY(2,IP),iz)
            v2  = BC3D_T(INDXY(3,IP),INDXY(2,IP),iz)
            v3  = BC3D_T(INDXY(1,IP),INDXY(4,IP),iz)
            v4  = BC3D_T(INDXY(3,IP),INDXY(4,IP),iz)
            ! if we have any fill values we have hit the bottom
            IF (u1.LT.FVP .OR. u2.LT.FVP .OR.&
                u3.LT.FVP .OR. u4.LT.FVP .OR.&
                v1.LT.FVP .OR. v2.LT.FVP .OR.&
                v3.LT.FVP .OR. v4.LT.FVP) THEN
               ADC3D_U(iz,IP) = FV
               ADC3D_V(iz,IP) = FV
            ELSE
               ! grab previous timestep for trapezoidal ule
               IF (iz.GT.1) THEN
                  UUAVG = UU
                  VVAVG = VV
               ENDIF
               UU = u1*WEIGHTS(1,IP) + u2*WEIGHTS(2,IP) +&
                    u3*WEIGHTS(3,IP) + u4*WEIGHTS(4,IP)
               VV = v1*WEIGHTS(1,IP) + v2*WEIGHTS(2,IP) +&
                    v3*WEIGHTS(3,IP) + v4*WEIGHTS(4,IP)
               ! save these values
               ADC3D_U(iz,IP) = UU
               ADC3D_V(iz,IP) = VV
               IF (iz.GT.1) THEN
                  UUAVG = 0.5d0*(UUAVG + UU)
                  VVAVG = 0.5d0*(VVAVG + VV)
                  dz = BC3D_Z(iz) - BC3D_Z(iz-1)
                  ADC2D_U(IP) = ADC2D_U(IP) + dz*UUAVG
                  ADC2D_V(IP) = ADC2D_V(IP) + dz*VVAVG
               ENDIF
               truebottom = iz
            ENDIF
         ENDDO
         ! depth average the velocities
         IF (truebottom.GT.1) THEN
            ADC2D_U(IP) = ADC2D_U(IP)/BC3D_Z(truebottom)
            ADC2D_V(IP) = ADC2D_V(IP)/BC3D_Z(truebottom)
         ELSE
            ADC2D_U(IP) = 0d0
            ADC2D_V(IP) = 0d0
         ENDIF
      ENDDO
      !
      ! loop through again to calculate Duu, Dvv, Duv
      !
      DO IP = 1,NP
         ! find predetermined max depth index
         idzmax = INDZ(IP)
         ! if we are too shallow skip
         IF ( idzmax.LT.2 ) CYCLE
         ! in case there is a bathy mismatch
         DO iz = 1,idzmax
            ! save previous vdiff, udiff for trapezoidal rule
            IF (iz.GT.1) THEN
               Udiffavg = Udiff
               Vdiffavg = Vdiff
            ENDIF
            ! if we hit the bottom exit this loop
            IF (ADC3D_U(iz,IP).LT.FVP .OR. ADC3D_V(iz,IP).LT.FVP) THEN
               EXIT
            ENDIF
            ! find Udiff, Vdiff
            Udiff = ADC3D_U(iz,IP) - ADC2D_U(IP)
            Vdiff = ADC3D_V(iz,IP) - ADC2D_V(IP)
            ! integrate
            IF (iz.GT.1) THEN
               dz = BC3D_Z(iz) - BC3D_Z(iz-1)
               Udiffavg = 0.5d0*(Udiff+Udiffavg)
               Vdiffavg = 0.5d0*(Vdiff+Vdiffavg)
               DUU_ADC(IP) = DUU_ADC(IP) + Udiffavg*Udiffavg*dz
               DVV_ADC(IP) = DVV_ADC(IP) + Vdiffavg*Vdiffavg*dz
               DUV_ADC(IP) = DUV_ADC(IP) + Udiffavg*Vdiffavg*dz
            ENDIF
         ENDDO
      ENDDO
      !
      ! Calculate momentum dispersion
      !
      DO IE = 1,NE
         NM1 = NM(IE,1)
         NM2 = NM(IE,2)
         NM3 = NM(IE,3)
         ddx1 = Dphi1Dx(IE)
         ddx2 = Dphi2Dx(IE)
         ddx3 = Dphi3Dx(IE)
         ddy1 = Dphi1Dy(IE)
         ddy2 = Dphi2Dy(IE)
         ddy3 = Dphi3Dy(IE)
         ! calculate this element's derivatives
         DUUx = DUU_ADC(NM1)*ddx1 +&
                DUU_ADC(NM2)*ddx2 +&
                DUU_ADC(NM3)*ddx3
         DVVy = DVV_ADC(NM1)*ddy1 +&
                DVV_ADC(NM2)*ddy2 +&
                DVV_ADC(NM3)*ddy3
         DUVx = DUV_ADC(NM1)*ddx1 +&
                DUV_ADC(NM2)*ddx2 +&
                DUV_ADC(NM3)*ddx3
         DUVy = DUV_ADC(NM1)*ddy1 +&
                DUV_ADC(NM2)*ddy2 +&
                DUV_ADC(NM3)*ddy3
        ! add this elements contribution to each node
        DispX_ADC(NM1) = DispX_ADC(NM1) + (DUUx+DUVy)*AREAS(IE)
        DispX_ADC(NM2) = DispX_ADC(NM2) + (DUUx+DUVy)*AREAS(IE)
        DispX_ADC(NM2) = DispX_ADC(NM2) + (DUUx+DUVy)*AREAS(IE)
        DispY_ADC(NM1) = DispY_ADC(NM1) + (DUVx+DVVy)*AREAS(IE)
        DispY_ADC(NM2) = DispY_ADC(NM2) + (DUVx+DVVy)*AREAS(IE)
        DispY_ADC(NM2) = DispY_ADC(NM2) + (DUVx+DVVy)*AREAS(IE)
      ENDDO
      !
      ! depth average
      !
      DO IP = 1,NP
         IF (DP(IP).GT.1d-3) THEN
            DispX_ADC(IP) = DispX_ADC(IP)/DP(IP)/TotalArea(IP)
            DispY_ADC(IP) = DispY_ADC(IP)/DP(IP)/TotalArea(IP)
         ELSE
            DispX_ADC(IP) = 0d0
            DispY_ADC(IP) = 0d0
         ENDIF
      ENDDO
      ! for memory management deallocate this stuff
      DEALLOCATE( ADC3D_U, ADC3D_V, ADC2D_U, ADC2D_V, &
                  DUU_ADC, DVV_ADC, DUV_ADC )
!-----------------------------------------------------------------------
      END SUBROUTINE Calc_BC2D_UV_ADCIRC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E  C A L C _ S T E R I C _ A D J U S T M E N T
!-----------------------------------------------------------------------
!     This subroutine takes the salinity and temperature read from GOFS
!     data and:
!        1) Converts from practical salinity to absolute salinity and
!           from temperature to conservative temperature on the GOFS
!           grid.
!        2) Calculates dynamic height anomaly w.r.t the bottom for each
!           ADCIRC node. 
!        3) Calculates baroclinic pressure anomaly:
!               P'(z) = P(z) - (1/D)\int_{-D}^{0}P(z)dz
!        4) Estimates the baroclinc sea level (BCSL_ADC) for each node
!-----------------------------------------------------------------------
      SUBROUTINE Calc_Steric_Adjustment()
      IMPLICIT NONE
      REAL*8 :: RHO(NZ), RHOAVG, DZ, DHA_avg, dha_step
      REAL*8 :: SA(NZ), CT(NZ), LAT(NZ), DHA(NZ)
      REAL(SZ) :: FV = -3D4, FVP = -3D4 + 1D-3, VS = 1D-3 
      ! persistent logical for if we are on the first call
      LOGICAL,SAVE :: FirstCall = .TRUE.
      ! looping variables for do loops
      INTEGER :: IP, IE, iz, idzmax, truebottom, NM1, NM2, NM3
      ! for calculating baroclinic pressures
      REAL*8 :: xx, yy, sp1, sp2, sp3, sp4, t1, t2, t3, t4,&
                sa1, sa2, sa3, sa4, ct1, ct2, ct3, ct4
      ! allocate output arrays if on the first call
      WRITE(6,*) "INFO: calculating BCSL"
      IF (FirstCall) THEN
         ALLOCATE( BCSL_ADC(NP) ) 
         FirstCall = .FALSE.
      ENDIF
      ! reset all values to 0
      BCSL_ADC(:)   = 0d0
      !
      ! get BCP_ADC, SigTS_ADC, NB_ADC, and NM_ADC
      !
      DO IP = 1,NP
         ! save lat and lon as yy,xx for readability
         xx = SLAM(IP)
         IF (minval(BC3D_Lon).GE.0d0.AND.xx.lt.0d0) xx = xx + 360d0

         yy = SFEA(IP)
         ! get max depth we are using
         idzmax = INDZ(IP)
         ! if we are only at the surface skip this whole loop
         IF (idzmax.lt.2) CYCLE
         ! reset the variable to track if we are at the bottom to 0
         truebottom = 0
         ! for my sanity reset the arrays for SA and CT to 0
         SA = 0d0
         CT = 0d0
         RHO = 0d0
         ! we also need a latitude array for the buoyancy frequency profile
         ! calculations
         Lat = yy
         DO iz = 1,idzmax
            ! for readability save the salinity and temps as individual values
            ! rather than using the indices
            ! salinity
            sp1 = BC3D_SP(INDXY(1,IP),INDXY(2,IP),iz)
            sp2 = BC3D_SP(INDXY(3,IP),INDXY(2,IP),iz)
            sp3 = BC3D_SP(INDXY(1,IP),INDXY(4,IP),iz)
            sp4 = BC3D_SP(INDXY(3,IP),INDXY(4,IP),iz)
            ! temperature
            t1  = BC3D_T(INDXY(1,IP),INDXY(2,IP),iz)
            t2  = BC3D_T(INDXY(3,IP),INDXY(2,IP),iz)
            t3  = BC3D_T(INDXY(1,IP),INDXY(4,IP),iz)
            t4  = BC3D_T(INDXY(3,IP),INDXY(4,IP),iz)
            ! if any of these are a fill value then we have hit the bottom and
            ! should use a fill value
            IF (sp1.lt.FVP .OR. sp2.lt.FVP .OR.&
                sp3.lt.FVP .OR. sp4.lt.FVP .OR.&
                t1.lt.FVP  .OR. t2.lt.FVP .OR.&
                t3.lt.FVP  .OR. t4.lt.FVP) THEN
            ELSE
               ! convert to absolute salinity
               sa1 = gsw_SA_from_SP(max(2d0,sp1),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(1,IP)),BC3D_LAT(INDXY(2,IP)))
               !
               sa2 = gsw_SA_from_SP(max(2d0,sp2),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(3,IP)),BC3D_LAT(INDXY(2,IP)))
               !
               sa3 = gsw_SA_from_SP(max(2d0,sp3),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(1,IP)),BC3D_LAT(INDXY(4,IP)))
               !
               sa4 = gsw_SA_from_SP(max(2d0,sp4),BC3D_Z(iz),&
                            BC3D_Lon(INDXY(3,IP)),BC3D_LAT(INDXY(4,IP)))
               ! convert to conservative temperature
               ct1 = gsw_CT_from_T(sa1,dble(t1),BC3D_Z(iz))
               ct2 = gsw_CT_from_T(sa2,dble(t2),BC3D_Z(iz))
               ct3 = gsw_CT_from_T(sa3,dble(t3),BC3D_Z(iz))
               ct4 = gsw_CT_from_T(sa4,dble(t4),BC3D_Z(iz))
               ! interpolate to our ADCIRC node
               SA(iz) = sa1*WEIGHTS(1,IP) + sa2*WEIGHTS(2,IP) +&
                        sa3*WEIGHTS(3,IP) + sa4*WEIGHTS(4,IP)
               CT(iz) = ct1*WEIGHTS(1,IP) + ct2*WEIGHTS(2,IP) +&
                        ct3*WEIGHTS(3,IP) + ct4*WEIGHTS(4,IP)
               RHO(iz) = gsw_rho(SA(iz),CT(iz),BC3D_Z(iz))
               truebottom = iz
            ENDIF
         ENDDO
         ! if applicable get the dynamic height anomaly
         IF (truebottom.GT.0) THEN
            DHA(:) = 0d0
            DHA(1:truebottom) = gsw_geo_strf_dyn_height(&
                                SA(1:truebottom), CT(1:truebottom),&
                                BC3D_Z(1:truebottom), &
                                BC3D_Z(truebottom))
            DHA_avg = 0d0
            DO iz = 1,truebottom
               IF (iz.EQ.1) THEN
                  dha_step = DHA(iz)
                  CYCLE
               ENDIF
               dha_step = (1d0/2d0)*(dha_step + (DHA(iz)))
               DZ = BC3D_Z(iz) - BC3D_Z(iz-1)
               DHA_avg = DHA_avg + dha_step*DZ
               dha_step = DHA(iz)
            ENDDO
            BCSL_adc(IP) = (DHA(1) - DHA_avg/BC3D_Z(truebottom))
            BCSL_adc(IP) = BCSL_adc(IP)/(g)
            IF (BCSL_adc(IP).GT.1d5) THEN
               BCSL_adc(IP) = 0
            ENDIF
         ENDIF
      ENDDO
!-----------------------------------------------------------------------
      END SUBROUTINE Calc_Steric_Adjustment
!-----------------------------------------------------------------------
      function haversine(deglon1,deglon2,deglat1,deglat2) result (dist)
          ! great circle distance -- adapted from Matlab 
          real*8,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real*8 :: a,c,dist,dlat,dlon,lat1,lat2
 
          dlat = deg2rad*(deglat2-deglat1)
          dlon = deg2rad*(deglon2-deglon1)
          lat1 = deg2rad*(deglat1)
          lat2 = deg2rad*(deglat2)
          a = ( sin(0.5d0*dlat) )**2 + &
                cos(lat1) * cos(lat2) * ( sin(0.5d0*dlon) )**2
          c = 2d0*asin( sqrt(a) )
          dist = Rearth*c
      end function haversine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T _ N E T C D F
!-----------------------------------------------------------------------
      subroutine initNetCDF()
      implicit none
      logical :: Fexist 
      integer :: ierr, data_dims(3), TEE 
      type(timedelta) :: DT
      type(datetime)  :: TNDT
      character(len=16) :: TSS

      if (myProc.eq.0) then   
         inquire(file=BC2D_Name,exist=Fexist)
      else
         Fexist = .true.
      endif
      if (.not.Fexist.and.OutType.LT.3) then
         ! Open file  
         call check_err(nf90_create(BC2D_Name, ncformat, ncid))
         ! Define dimensions
         call check_err(nf90_def_dim(ncid, 'time', nf90_unlimited, &
                        timenc_dim_id))
         call check_err(nf90_def_dim(ncid, 'strlen', 16, strlen_dim_id))
         call check_err(nf90_def_dim(ncid,'NX',NX,NX_dim_id))
         call check_err(nf90_def_dim(ncid,'NY',NY,NY_dim_id))
         if (OutType.gt.0) &
           call check_err(nf90_def_dim(ncid,'NYY',NY-1,NYY_dim_id))
         ! CPB 01/21/2022: changed condition to always define NYYY_dim_id because of netcdf error
         ! This should not be an issue since it just results in the latitude for Disp to be written
         ! even if it is not used in the adcirc run
         if (OutType.gt.0) &
           call check_err(nf90_def_dim(ncid,'NYYY',NY-2,NYYY_dim_id))
         ! Define vars 
         data_dims = [strlen_dim_id, timenc_dim_id, 1]
         call def_var_att(ncid,'time',nf90_char, data_dims(1:2), &
                          timenc_id,'UTC datetime','YYYY-MM-DD HH:mm')
         call def_var_att(ncid,'lon',nf90_double, [NX_dim_id], &
                          lon_id,'longitude','degrees')
         call def_var_att(ncid,'lat',nf90_double, [NY_dim_id], &
                          lat_id,'latitude','degrees')
         if (OutType.gt.0) then
           call def_var_att(ncid,'lonc',nf90_double, [NX_dim_id], &
                          lonc_id,'longitude for BPGX','degrees')
           call def_var_att(ncid,'latc',nf90_double, [NYY_dim_id], &
                          latc_id,'latitude for BPGY','degrees')
           call def_var_att(ncid,'lats',nf90_double, [NYYY_dim_id], &
                          lats_id,'latitude for Disp','degrees')
           data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'BPGX',nf90_type, data_dims, BPGX_id, &
                        'east-west depth-averaged baroclinic pressure '&
                        //'gradient','ms^-2')
           data_dims = [NX_dim_id, NYY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'BPGY',nf90_type, data_dims, BPGY_id, &
                      'north-south depth-averaged baroclinic pressure '&
                        //'gradient','ms^-2')
           data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'SigTS',nf90_type, data_dims,   &
                          SigTS_id, 'surface sigmat density',  &
                          'kgm^-3',SigT0)
           call def_var_att(ncid,'MLD',nf90_type, data_dims, MLD_id, &
                          'mixed-layer depth ratio','[]',DFV)
         !                 'mixed-layer depth','m',DFV)
         endif
         data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
         call def_var_att(ncid,'NB',nf90_type, data_dims, NB_id, &
                         'buoyancy frequency at the seabed','s^-1')
         call def_var_att(ncid,'NM',nf90_type, data_dims, NM_id, &
                         'depth-averaged buoyancy frequency','s^-1')
         if (OutType.eq.2) then
            call def_var_att(ncid,'KE',nf90_type, data_dims, KE_id, &
                         'depth-averaged kinetic energy','m^2s^-2')
            data_dims = [NX_dim_id, NYYY_dim_id, timenc_dim_id]
            call def_var_att(ncid,'CDisp',nf90_type, data_dims, CD_id, &
                   'momentum dispersion quadratic bottom friction','[]')
            call def_var_att(ncid,'DispX',nf90_type, data_dims, &
                  DispX_id, 'depth-averaged x-momentum dispersion',&
                  'ms^-2')
            call def_var_att(ncid,'DispY',nf90_type, data_dims, &
                  DispY_id, 'depth-averaged y-momentum dispersion',&
                  'ms^-2')
         endif
         ! Allowing vars to deflate
         call check_err(nf90_def_var_deflate(ncid, lon_id, 1, 1, dfl))
         call check_err(nf90_def_var_deflate(ncid, lat_id, 1, 1, dfl))
         if (OutType.gt.0) then
            call check_err(nf90_def_var_deflate(ncid, lonc_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, latc_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, lats_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, BPGX_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, BPGY_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, SigTS_id,1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, MLD_id,  1,1,dfl))
         endif
         call check_err(nf90_def_var_deflate(ncid, NB_id, 1, 1, dfl))
         call check_err(nf90_def_var_deflate(ncid, NM_id, 1, 1, dfl))
         if (OutType.eq.2) then
            call check_err(nf90_def_var_deflate(ncid, KE_id,   1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, CD_id,   1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, DispX_id,1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, DispY_id,1,1,dfl))
         endif
         ! Close define mode
         call check_err(nf90_close(ncid))
         
         ! Put X, Y, Z on
         call check_err(nf90_open(BC2D_Name, nf90_write, ncid))
         if (is == 1) then
            ! Already -180/180 orientation
            call check_err(nf90_put_var(ncid, lon_id, BC3D_Lon))
            if (OutType.gt.0) &
               call check_err(nf90_put_var(ncid, lonc_id,   &
                 [0.5d0*(BC3D_Lon(1:NX-1) + BC3D_Lon(2:NX)),&
                  0.5d0*(BC3D_Lon(NX) + BC3D_Lon(1)+360d0)]))
         else
            ! Change to -180/180 orientation
            call check_err(nf90_put_var(ncid, lon_id,      &
                 [BC3D_Lon(is:NX)-360d0, BC3D_Lon(1:is-1)]))
            if (OutType.gt.0) &
               call check_err(nf90_put_var(ncid, lonc_id,             &
                 [0.5d0*(BC3D_Lon(is:NX-1) + BC3D_Lon(is+1:NX))-360d0,&
                  0.5d0*(BC3D_Lon(NX)-360d0 + BC3D_Lon(1)),           &
                  0.5d0*(BC3D_Lon(1:is-1)  + BC3D_Lon(2:is))]))
         endif
         call check_err(nf90_put_var(ncid, lat_id, BC3D_Lat))
         if (OutType.gt.0) &
            call check_err(nf90_put_var(ncid, latc_id, &
              0.5d0*(BC3D_Lat(1:NY-1) + BC3D_Lat(2:NY))))
         if (OutType.eq.2) &
            call check_err(nf90_put_var(ncid, lats_id,BC3D_Lat(2:NY-1)))
         call check_err(nf90_close(ncid))
      ELSEIF ( .not.Fexist.AND.OutType.GE.3 ) THEN
         CALL initNetCDF_adc()
      endif
      ! Barrier to ensure wait until netcdf is created by first
      ! processor
      call MPI_Barrier(MPI_Comm_World,ierr)     
      
      if (Fexist) then 
         ! Get the dim and var ids
         call check_err(nf90_open(BC2D_Name, nf90_nowrite, ncid))
         call check_err(nf90_inq_dimid(ncid,'time', timenc_dim_id))
         call Check_err(nf90_inquire_dimension(ncid,&
                        timenc_dim_id,len=TEE))
         call check_err(nf90_inq_varid(ncid,'time', timenc_id))
         ! Let's get the start time index
         do tsind = 1,TEE
            call check_err(nf90_get_var(ncid, timenc_id, &
                           TSS,start=[1, tsind],count=[16, 1]))
            TNDT = strptime(trim(TSS),"%Y-%m-%d %H:%M")
            DT   = TSDT - TNDT
            if (DT%total_seconds() <= 0) exit
         enddo
         if (myProc.eq.0) write(6,*) 'tsind = ',tsind
         if (OutType.gt.0) then
            call check_err(nf90_inq_varid(ncid,'BPGX', BPGX_id))
            call check_err(nf90_inq_varid(ncid,'BPGY', BPGY_id))
            call check_err(nf90_inq_varid(ncid,'SigTS', SigTS_id))
            call check_err(nf90_inq_varid(ncid,'MLD', MLD_id))
         endif
         call check_err(nf90_inq_varid(ncid,'NB', NB_id))    
         call check_err(nf90_inq_varid(ncid,'NM', NM_id))
         if (OutType.eq.2.OR.OutType.EQ.4) then
            !call check_err(nf90_inq_varid(ncid,'KE', KE_id))
            !call check_err(nf90_inq_varid(ncid,'CDisp', CD_id))
            call check_err(nf90_inq_varid(ncid,'DispX', DispX_id))
            call check_err(nf90_inq_varid(ncid,'DispY', DispY_id))
         endif
         call check_err(nf90_close(ncid))
      endif
      end subroutine initNetCDF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T N E T C D F _ ADC
!-----------------------------------------------------------------------
      subroutine initNetCDF_adc()
      implicit none
      logical :: Fexist 
      integer :: ierr, data_dims(3), TEE 
      type(timedelta) :: DT
      type(datetime)  :: TNDT
      character(len=16) :: TSS
      !
      ! print to screen for intormation
      !
      WRITE(6,*) 'INFO: initializing NetCDF for ADCIRC output'
      !
      ! Open file  
      !
      call check_err(nf90_create(BC2D_Name, ncformat, ncid))
      !
      ! Define dimensions
      !
      call check_err(nf90_def_dim(ncid, 'time', nf90_unlimited, &
                     timenc_dim_id))
      call check_err(nf90_def_dim(ncid, 'strlen', 16, strlen_dim_id))
      call check_err(nf90_def_dim(ncid,'node',NP,node_dim_id))
      call check_err(nf90_def_dim(ncid,'nele',NE,nele_dim_id))
      call check_err(nf90_def_dim(ncid,'nvertex',3,nvertex_dim_id))
      !
      ! define variables
      !
      data_dims = [strlen_dim_id, timenc_dim_id, 1]
      !
      ! define time variable
      !
      call def_var_att(ncid,'time',nf90_char, data_dims(1:2), &
                       timenc_id,'UTC datetime','YYYY-MM-DD HH:mm')
      !
      ! define mesh variables (x, y, element, b)
      !
      call def_var_att(ncid,'x',nf90_double, [node_dim_id], &
                       lon_id,'longitude','degrees')
      call def_var_att(ncid,'y',nf90_double, [node_dim_id], &
                          lat_id,'latitude','degrees')
      call def_var_att(ncid,'depth',nf90_double, [node_dim_id], &
                          depth_id,'distance below geoid','m')
      call def_var_att(ncid,'element',nf90_int,&
                          [nele_dim_id, nvertex_dim_id], &
                          ele_id,'element','nondimensional')
      !
      ! define calculated values
      !
      data_dims = [node_dim_id,timenc_dim_id,1]
      call def_var_att(ncid,'BPGX',nf90_type, data_dims(1:2), BPGX_id, &
                     'east-west depth-averaged baroclinic pressure '&
                     //'gradient','ms^-2')
      call def_var_att(ncid,'BPGY',nf90_type, data_dims(1:2), BPGY_id, &
                 'north-south depth-averaged baroclinic pressure '&
                   //'gradient','ms^-2')
      call def_var_att(ncid,'SigTS',nf90_type, data_dims(1:2),   &
                     SigTS_id, 'surface sigmat density',  &
                     'kgm^-3',SigT0)
      call def_var_att(ncid,'MLD',nf90_type, data_dims(1:2), MLD_id, &
                     'mixed-layer depth ratio','[]',DFV)
      call def_var_att(ncid,'NB',nf90_type, data_dims(1:2), NB_id, &
                      'buoyancy frequency at the seabed','s^-1')
      call def_var_att(ncid,'NM',nf90_type, data_dims(1:2), NM_id, &
                      'depth-averaged buoyancy frequency','s^-1')
      !
      ! define momentum dispersion if appropriate
      !
      IF (OutType.GT.3) THEN
         call def_var_att(ncid,'DispX',nf90_type, data_dims(1:2),&
                          DispX_id, &
                         'depth-averaged x-momentum dispersion','ms^-2')
         call def_var_att(ncid,'DispY',nf90_type, data_dims(1:2),&
                          DispY_id, &
                         'depth-averaged x-momentum dispersion','ms^-2')
      ENDIF
      !
      ! Allowing vars to deflate
      !
      call check_err(nf90_def_var_deflate(ncid, lon_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, lat_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, depth_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, ele_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, BPGX_id, 1,1,dfl))
      call check_err(nf90_def_var_deflate(ncid, BPGY_id, 1,1,dfl))
      call check_err(nf90_def_var_deflate(ncid, SigTS_id,1,1,dfl))
      call check_err(nf90_def_var_deflate(ncid, MLD_id,  1,1,dfl))
      call check_err(nf90_def_var_deflate(ncid, NB_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, NM_id, 1, 1, dfl))
      IF (OutType.GT.3) THEN
         call check_err(nf90_def_var_deflate(ncid, DispX_id, 1, 1, dfl))
         call check_err(nf90_def_var_deflate(ncid, DispY_id, 1, 1, dfl))
      ENDIF
      ! 
      ! put on mesh variables
      !
      call check_err(nf90_put_var(ncid, lon_id, SLAM))
      call check_err(nf90_put_var(ncid, lat_id, SFEA))
      call check_err(nf90_put_var(ncid, depth_id, DP))
      call check_err(nf90_put_var(ncid, ele_id, NM))
      !
      ! close file
      !
      call check_err(nf90_close(ncid))
      end subroutine initNetCDF_adc
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    D E F _ V A R _ A T T
!-----------------------------------------------------------------------
      subroutine def_var_att(ncid, Sname, NFtype, dims, var_id, Lname, &
                 Units, FillValue, Factor, Offset)
      implicit none
      integer,intent(in)  :: ncid, NFtype 
      integer,intent(in),dimension(:) :: dims
      integer,intent(out) :: var_id   
      character(*),intent(in) :: Sname, Lname, Units 
      real(sz),intent(in),optional :: FillValue, Offset, Factor

      call check_err(nf90_def_var(ncid, Sname, NFtype, dims, var_id))
      call check_err(nf90_put_att(ncid, var_id, 'long_name', Lname))
      call check_err(nf90_put_att(ncid, var_id, 'units', Units))
      
      if (present(FillValue)) then
         call Check_err(nf90_put_att(ncid, var_ID, &
                        '_FillValue',FillValue))
        !call Check_err(nf90_put_att(ncid, var_ID,'scale_factor',Factor))
        !call Check_err(nf90_put_att(ncid, var_ID,'add_offset',Offset))
      endif
!
      end subroutine def_var_att       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    U P D A T E N E T C D F
!-----------------------------------------------------------------------
      subroutine UpdateNetCDF(IT,OKflag)
      implicit none
      integer, intent(in) :: IT
      logical, intent(in) :: OKflag
      integer,dimension(3) :: start, kount, kountX, kountY, kountYY
      integer  :: startT, mpistat(mpi_status_size)
      IF (OutType.EQ.5) THEN
         CALL writeBCSLNC(IT,OKflag)
         return
      ENDIF 
      ! Make netcdf/get dimension and variable IDs 
      write(6,*) 'IT = ',IT
      if (IT.eq.0) then 
         call initNetCDF()
         write(6,*)'PROC ',myproc,'initiated netcdf'
      endif
      !start = [1, 1, mnProc*IT + myProc + tsind];
      IF (OutType.LT.3) THEN
         kount = [NX, NY, 1];
         kountY = [NX, NY-1, 1];
         kountYY = [NX, NY-2, 1];
      ELSE
         kount = [1,NP,1]
      ENDIF
      if ( (IT > 0 .or. myProc > 0) .and. mnProc > 1) then
         ! Receiving message from prev processor to start writing
         call mpi_recv(startT,1,mpi_integer,prevProc,0,MPI_Comm_World,&
                       mpistat,ierr)
      endif

      if (OKflag) then
         if ((IT.eq.0.and.myProc.eq.0).or.mnProc.eq.1) then
            startT = tsind + IT
         else
            startT = startT + 1
         endif
         start = [1, 1, startT];
         write(6,*) 'UpdateNetCDF:',myProc,&
                    CurDT%strftime("%Y-%m-%d %H:%M")
         call check_err(nf90_open(BC2D_Name, nf90_write, ncid))
         call check_err(nf90_put_var(ncid, timenc_id,    &
                        CurDT%strftime("%Y-%m-%d %H:%M"),&
                        [1, start(3)],[16, kount(3)]))
         IF (OutType.LT.3) THEN
            if (OutType.gt.0) then
               call put_var(ncid, BPGX_id,BC2D_BX, start, kount)
               call put_var(ncid, BPGY_id, BC2D_BY, start, kountY)
               call put_var(ncid, SigTS_id, BC2D_SigTS, start, kount)
               call put_var(ncid, MLD_id, BC2D_MLD, start, kount)
            endif
            call put_var(ncid, NB_id, BC2D_NB, start, kount)
            call put_var(ncid, NM_id, BC2D_NM, start, kount)
            if (OutType.eq.2) then
               call put_var(ncid, KE_id, BC2D_KE, start, kount)
               call put_var(ncid, CD_id, BC2D_CD, start, kountYY)
               call put_var(ncid, DispX_id, BC2D_DispX, start, kountYY)
               call put_var(ncid, DispY_id, BC2D_DispY, start, kountYY)
            endif
         ELSEIF (OutType.GE.3) THEN
            call put_var_adc(ncid, BPGX_id, BPG_ADCx, start(2:3), &
                         kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied BPGX'
            call put_var_adc(ncid, BPGY_id, BPG_ADCy, start(2:3), & 
                         kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied BPGY'
            call put_var_adc(ncid, SigTS_id, SigTS_ADC, start(2:3),&
                         kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied SigTS'
            call put_var_adc(ncid, MLD_id, MLD_ADC, start(2:3), &
                             kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied MLD'
            call put_var_adc(ncid, NB_id, NB_ADC, start(2:3), &
                             kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied NB'
            call put_var_adc(ncid, NM_id, NM_ADC, start(2:3), &
                             kount(2:3))
            write(6,*) 'PROC',myProc,' succesfully applied NM'
            IF (OutType.GT.3) THEN
               call put_var_adc(ncid, DispX_id, DispX_ADC, start(2:3), &
                                kount(2:3))
               call put_var_adc(ncid, DispY_id, DispY_ADC, start(2:3), &
                                kount(2:3))

            ENDIF
         ENDIF
         call check_err(nf90_close(ncid))
      else
         write(6,*) 'Skipping Update of NetCDF:',&
                     myProc,CurDT%strftime("%Y-%m-%d %H:%M")
      endif

      if (mnProc > 1) then
         ! Telling next processor to start writing 
         call mpi_send(startT,1,mpi_integer,nextProc,0,&
                       MPI_Comm_World,ierr)
      endif 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine UpdateNetCDF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    W R I T E B C S L N C
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine writeBCSLNC(IT,OKflag)
      implicit none
      integer, intent(in) :: IT
      logical, intent(in) :: OKflag
      integer,dimension(3) :: start, kount, kountX, kountY, kountYY
      integer  :: startT, mpistat(mpi_status_size)
      ! Make netcdf/get dimension and variable IDs 
      write(6,*) 'IT = ',IT
      if (IT.eq.0) then 
         call initNetCDF_BCSL()
         write(6,*)'PROC ',myproc,'initiated netcdf'
      endif
      kount = [1,NP,1]
      if ( (IT > 0 .or. myProc > 0) .and. mnProc > 1) then
         ! Receiving message from prev processor to start writing
         call mpi_recv(startT,1,mpi_integer,prevProc,0,MPI_Comm_World,&
                       mpistat,ierr)
      endif

      if (OKflag) then
         if ((IT.eq.0.and.myProc.eq.0).or.mnProc.eq.1) then
            startT = tsind + IT
         else
            startT = startT + 1
         endif
         start = [1, 1, startT];
         write(6,*) 'UpdateNetCDF:',myProc,&
                    CurDT%strftime("%Y-%m-%d %H:%M")
         call check_err(nf90_open(BC2D_Name, nf90_write, ncid))
         call check_err(nf90_put_var(ncid, timenc_id,    &
                        CurDT%strftime("%Y-%m-%d %H:%M"),&
                        [1, start(3)],[16, kount(3)]))
         call put_var_adc(ncid, BCSL_id, BCSL_ADC, start(2:3), &
                      kount(2:3))
         call check_err(nf90_close(ncid))
      else
         write(6,*) 'Skipping Update of NetCDF:',&
                     myProc,CurDT%strftime("%Y-%m-%d %H:%M")
      endif

      if (mnProc > 1) then
         ! Telling next processor to start writing 
         call mpi_send(startT,1,mpi_integer,nextProc,0,&
                       MPI_Comm_World,ierr)
      endif 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine writeBCSLNC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T N E T C D F _ B C S L
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine initNetCDF_BCSL()
      implicit none
      logical :: Fexist 
      integer :: ierr, data_dims(3), TEE 
      type(timedelta) :: DT
      type(datetime)  :: TNDT
      character(len=16) :: TSS
      !
      ! print to screen for intormation
      !
      WRITE(6,*) 'INFO: initializing NetCDF for Baroclinic Sea Level'
      !
      ! Open file  
      !
      call check_err(nf90_create(BC2D_Name, ncformat, ncid))
      !
      ! Define dimensions
      !
      call check_err(nf90_def_dim(ncid, 'time', nf90_unlimited, &
                     timenc_dim_id))
      call check_err(nf90_def_dim(ncid, 'strlen', 16, strlen_dim_id))
      call check_err(nf90_def_dim(ncid,'node',NP,node_dim_id))
      call check_err(nf90_def_dim(ncid,'nele',NE,nele_dim_id))
      call check_err(nf90_def_dim(ncid,'nvertex',3,nvertex_dim_id))
      !
      ! define variables
      !
      data_dims = [strlen_dim_id, timenc_dim_id, 1]
      !
      ! define time variable
      !
      call def_var_att(ncid,'time',nf90_char, data_dims(1:2), &
                       timenc_id,'UTC datetime','YYYY-MM-DD HH:mm')
      !
      ! define mesh variables (x, y, element, b)
      !
      call def_var_att(ncid,'x',nf90_double, [node_dim_id], &
                       lon_id,'longitude','degrees')
      call def_var_att(ncid,'y',nf90_double, [node_dim_id], &
                          lat_id,'latitude','degrees')
      call def_var_att(ncid,'depth',nf90_double, [node_dim_id], &
                          depth_id,'distance below geoid','m')
      call def_var_att(ncid,'element',nf90_int,&
                          [nele_dim_id, nvertex_dim_id], &
                          ele_id,'element','nondimensional')
      !
      ! define calculated values
      !
      data_dims = [node_dim_id,timenc_dim_id,1]
      call def_var_att(ncid,'BCSL',nf90_type, data_dims(1:2), BCSL_id, &
                     'baroclinic sea level ','m')
      !
      ! Allowing vars to deflate
      !
      call check_err(nf90_def_var_deflate(ncid, lon_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, lat_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, depth_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, ele_id, 1, 1, dfl))
      call check_err(nf90_def_var_deflate(ncid, BCSL_id, 1,1,dfl))
      ! 
      ! put on mesh variables
      !
      call check_err(nf90_put_var(ncid, lon_id, SLAM))
      call check_err(nf90_put_var(ncid, lat_id, SFEA))
      call check_err(nf90_put_var(ncid, depth_id, DP))
      call check_err(nf90_put_var(ncid, ele_id, NM))
      !
      ! close file
      !
      call check_err(nf90_close(ncid))
      end subroutine initNetCDF_BCSL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    P U T _ V A R
!-----------------------------------------------------------------------
      subroutine put_var(ncid, var_id, InArray, start, kount)
      implicit none
      integer,intent(in)  :: ncid, var_id, start(3), kount(3) 
      real(sz),intent(in) :: InArray(:,:)
      real(sz),allocatable :: Temp2D(:,:)

      if (is > 1) then
         ! If 0/360, re-orientate
         allocate(Temp2D(kount(1),kount(2)))
         Temp2D(1:kount(1)-is+1,:) = InArray(is:kount(1),:)
         Temp2D(kount(1)-is:kount(1),:)  = InArray(1:is-1,:)
         ! Write out into netcdf
         call check_err(nf90_put_var(ncid, var_id, Temp2D,  &
                        start, kount))
      else
         ! Write outinto netcdf
         call check_err(nf90_put_var(ncid, var_id, InArray, &
                        start, kount))
      endif
!
      end subroutine put_var       
!-----------------------------------------------------------------------
!     S U B R O U T I N E    P U T _ V A R _ ADC
!-----------------------------------------------------------------------
      subroutine put_var_adc(ncid, var_id, InArray, start, kount)
      implicit none
      integer,intent(in)  :: ncid, var_id, start(2), kount(2) 
      real(sz),intent(in) :: InArray(:)
      ! Write outinto netcdf
      call check_err(nf90_put_var(ncid, var_id, InArray, &
                     start, kount))
!
      end subroutine put_var_adc
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ F 1 4
!-----------------------------------------------------------------------
!     CPB: reads in lon, lat, depth, and connectivity table from ADCIRC
!     grid. Does not read in any boundary information.
!-----------------------------------------------------------------------
      subroutine Read_F14()
      implicit none
      logical :: Fexist
      INTEGER :: errorio
      ! looping variables
      INTEGER :: ii, dmy, dmy2
!
      ! open fort.14 file
      OPEN(14,FILE=trim(fort14),STATUS='OLD',ACTION='READ',&
      IOSTAT=errorio)
      IF (errorio.ne.0) THEN
         write(6,*) 'fort.14 file could not be opened' 
         call MPI_Abort(MPI_Comm_World,0,ierr)
      ENDIF
      READ(14,*)        ! skip First line
      READ(14,*) NE, NP
      ! allocate things
      ALLOCATE(SLAM(NP), SFEA(NP), DP(NP), NM(NE,3))
      ! read in lat, lon, depth
      DO ii = 1,NP
         READ(14,*) dmy, slam(ii), sfea(ii), dp(ii)
      ENDDO
      ! read in connectivity table
      DO ii = 1,NE
         read(14,*) dmy, dmy2, NM(ii,1), NM(ii,2), NM(ii,3)
      ENDDO
      CLOSE(14)
      WRITE(6,*) 'INFO: Successfully read fort.14 file'
!
      end subroutine Read_F14
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE Calc_Areas()
!-----------------------------------------------------------------------
!     This calculates the areas of each element (Areas) as well as the
!     total area connected to each node (TotalArea). This is slightly
!     different than ADCIRC because ADCIRC saves 2*area in Areas
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! dummy variables for holding elemental nodes
      REAL*8,DIMENSION(3) :: LATVE, LONVE, XVE, YVE
      ! looping variable
      INTEGER :: ii, NMTEMP(3)
      ! Intermediate variables to make this more readable while building
      ! areas, etc.
      REAL*8 :: x1, x2, x3, y1, y2, y3, &
                x2mx1, x3mx2, x1mx3, &
                y2my1, y3my2, y1my3
      ! jacobian to check if ccw
      REAL*8 :: JAC1, DLED1

      ! allocate Areas, TotalArea, FDXE, FDYE 
      ALLOCATE( Areas(NE), TotalArea(NP), FDXE(3,NE), FDYE(3,NE) )
      TotalArea = 0d0
      DO ii = 1,NE
         ! lon/lat in degrees
         LONVE(1) = SLAM( NM(ii,1) )
         LONVE(2) = SLAM( NM(ii,2) )
         LONVE(3) = SLAM( NM(ii,3) )
         LATVE(1) = SFEA( NM(ii,1) )
         LATVE(2) = SFEA( NM(ii,2) )
         LATVE(3) = SFEA( NM(ii,3) )
         LONVE = MODULO(LONVE,360d0)
         CALL Cal_Jac(JAC1,LONVE,LATVE)
         CALL Cal_Edgelength( DLED1, LONVE, LATVE )
         IF ( JAC1.LT.0d0.AND.DLED1.LT.360) THEN
             NMTEMP = NM(ii,:)
             NM(ii,1) = NMTEMP(3)
             NM(ii,2) = NMTEMP(2)
             NM(ii,3) = NMTEMP(1)
         ENDIF
         ! lon/lat in degrees
         LONVE(1) = SLAM( NM(ii,1) )
         LONVE(2) = SLAM( NM(ii,2) )
         LONVE(3) = SLAM( NM(ii,3) )
         LATVE(1) = SFEA( NM(ii,1) )
         LATVE(2) = SFEA( NM(ii,2) )
         LATVE(3) = SFEA( NM(ii,3) )
         ! get projected coordinates
         CALL CAL_ELXV_SPCOOR(XVE, YVE, LONVE, LATVE)
         ! GET X1-X2, ETC
         x1 = XVE(1)
         x2 = XVE(2)
         x3 = XVE(3)
         y1 = YVE(1)
         y2 = YVE(2)
         y3 = YVE(3)
         ! edgelengths
         x2mx1 = x2 - x1
         x3mx2 = x3 - x2
         x1mx3 = x1 - x3
         y2my1 = y2 - y1
         y3my2 = y3 - y2
         y1my3 = y1 - y3
         ! store these for later use in calculating derivatives
         FDXE(1,ii) = -Y3mY2 
         FDXE(2,ii) = -Y1mY3 
         FDXE(3,ii) = -Y2mY1 
         FDYE(1,ii) =  X3mX2 
         FDYE(2,ii) =  X1mX3 
         FDYE(3,ii) =  X2mX1 
         ! store areas
         AREAS(ii)=0.5d0*((X1mX3)*(-Y3mY2)+(X3mX2)*(Y1mY3))
         ! get totalarea
         TotalArea(NM(ii,1)) = TotalArea(NM(ii,1)) + Areas(ii)
         TotalArea(NM(ii,2)) = TotalArea(NM(ii,2)) + Areas(ii)
         TotalArea(NM(ii,3)) = TotalArea(NM(ii,3)) + Areas(ii)
      ENDDO
      WRITE(6,*) 'INFO: successfully calculated areas'
!
      END SUBROUTINE Calc_Areas
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE Calc_Derivatives()
!-----------------------------------------------------------------------
!     This calculates the derivatives of the basis functions of each
!     element
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! dummy variables for holding elemental nodes
      REAL*8,DIMENSION(3) :: LATVE, LONVE, XVE, YVE
      ! looping variable
      INTEGER :: ii, NM1, NM2, NM3
      ! Intermediate variables to make this more readable while building
      ! areas, etc.
      REAL*8 :: sfdxfac, sfdyfac, AreaIE2, AreaEle
      ! allocate memory for derivatives and adjustment factors for
      ! derivative calculations
      ALLOCATE( Dphi1Dx(NE), Dphi2Dx(NE), Dphi3Dx(NE), &
                Dphi1Dy(NE), Dphi2Dy(NE), Dphi3Dy(NE) )
      ! these are used in the calculation of the cylindrical projection
      ! derivatives and will be deallocated later this subroutine
      ALLOCATE( SFMXEle(NE), SFMyEle(NE) )
      ! convert to radians for this
      SLAM = SLAM*DEG2RAD
      SFEA = SFEA*DEG2RAD
      !WRITE(6,*) 'about to enter Compute_Cylinproj'
      CALL COMPUTE_CYLINPROJ_SFAC( SLAM, SFEA )
      ! back to degrees
      SLAM = SLAM*RAD2DEG
      SFEA = SFEA*RAD2DEG
      ! now that we have adjustment factors calculate derivatives
      DO ii = 1,NE
         NM1 = NM(ii,1)
         NM2 = NM(ii,2)
         NM3 = NM(ii,3)
         sfdxfac = SFmxEle(ii)
         sfdyfac = SFmyEle(ii)
         AreaEle = Areas(ii)
         AreaIE2 = 2*AreaEle
         ! calculate the derivatives
         Dphi1Dx(ii) = FDXE(1,ii)*sfdxfac/AreaIE2
         Dphi2Dx(ii) = FDXE(2,ii)*sfdxfac/AreaIE2
         Dphi3Dx(ii) = FDXE(3,ii)*sfdxfac/AreaIE2
         Dphi1Dy(ii) = FDYE(1,ii)*sfdyfac/AreaIE2
         Dphi2Dy(ii) = FDYE(2,ii)*sfdyfac/AreaIE2
         Dphi3Dy(ii) = FDYE(3,ii)*sfdyfac/AreaIE2
      ENDDO
      ! deallocate thngs we don't need anymore
      DEALLOCATE( SFmxEle, SFmyEle, FDXE, FDYE )
      WRITE(6,*) 'INFO: Found Dphi_i/Dx_i'
!
      END SUBROUTINE Calc_Derivatives
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE CAL_ELXV_SPCOOR( XVE, YVE, LONVE, LATVE )
!-----------------------------------------------------------------------
!     This is pulled directly from the ADCIRC source code with some
!     modifications
!    
!     INPUTS:
!      - LONVE = 3 longitude values of the nodes of an element in
!                degrees
!      - LATVE = 3 latitude values of the nodes of an element in degrees
!
!     OUTPUTS:
!      - XVE = 3 mapped coordinates of the nodes of an element
!      - YVE = 3 mapped coordinates of the nodes of an element
!
!     NOTES: I modified this to always use Mercator rather than have
!     options since it should be fine for this pre-processor. I also
!     changed it so that you don't have to swap back and forth between
!     degrees and radians quite so much.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:):: XVE, YVE, LONVE, LATVE
      INTEGER :: ELEID 
      !c Local 
      INTEGER:: II, IDX, SPF
      REAL*8:: XC, YC
      REAL*8:: Jac1, Jac2
      REAL*8, dimension(3):: LONM, LATM, LONTMP, LATTMP
      REAL*8, dimension(3):: DLX, DL
      REAL*8 :: DLED1, DLED2
  
      LATM(:) = LATVE(1:3)
      LATTMP(:) = LATVE(1:3)
      !c Adjust so that 0 <= lon <= 360 
      LONM(:) = MODULO( LONVE(1:3), 360.0d0 ) ; 
      LONTMP = LONM ; 
      CALL CAL_JAC( Jac1, LONM, LATM )
      CALL CAL_EDGELENGTH( DLED1, LONM, LATM ) ; 
   
      SPF = 0 ; 
      IF (  Jac1 < 0.0d0 .OR. &
          (Jac1 > 0.0d0 .AND. DLED1 > 360.0d0 ) ) THEN
         !c An element has a circular index. 
         !c Put a wrapped node on the other side
         !c A wrapped node is on the left of 180
         IF ( COUNT( LONM > 180.0d0 ) == 1 ) THEN
            IDX = sum( merge( (/ 1, 2, 3 /), &
                (/ 0, 0, 0 /), LONM > 180.0d0 ) ) ;
            LONM(IDX) = LONM(IDX) - 360.0d0 ; 
         END IF
         
         !c A wrapped node is on the right of 180
         IF ( COUNT( LONM < 180.0d0 ) == 1 ) THEN
            IDX = sum( merge( (/ 1, 2, 3 /), &
                (/ 0, 0, 0 /), LONM < 180.0d0 ) ) ;
            
            LONM(IDX) = LONM(IDX) + 360.0d0 ; 
         END IF        
         CALL CAL_JAC( Jac2, LONM, LATM )
         CALL CAL_EDGELENGTH( DLED2, LONM, LATM ) ; 
         SPF = 1 ;
         IF ( Jac2 < 0.0d0 .OR. &
             (Jac1 > 0.0d0 .AND. DLED2 > DLED1) ) THEN
            LONM = LONTMP ;  
            SPF = 0 ; 
         END IF
      END IF
            
      !c deg -> rad
      LONM = LONM*DEG2RAD ;
      LATM = LATM*DEG2RAD ;
      DO II = 1, 3
         ! set slam0 and sfea0 to just be 0 and 45 since that is what I
         ! do in my fort.15
         CALL CYLINDERMAP(XC, YC, LONM(II), LATM(II), slam0,sfea0)
         XVE(II) = XC ;
         YVE(II) = YC ; 
      END DO      
      
      RETURN ;
      END SUBROUTINE CAL_ELXV_SPCOOR
      
      SUBROUTINE CAL_EDGELENGTH( DEDGE, LON, LAT ) 
        IMPLICIT NONE
        REAL*8:: DEDGE, LON(3), LAT(3)
        !c local c!
        REAL*8, dimension(3):: DLX, DLY
        DLX = (/ LON(2) - LON(1), &
            LON(3) - LON(2), LON(1) - LON(3) /) ; 
        DLY = (/ LAT(2) - LAT(1), &
            LAT(3) - LAT(2), LAT(1) - LAT(3) /) ; 
        
        DEDGE = SQRT( SUM(DLX*DLX + DLY*DLY) ) ; 
        RETURN ;
      END SUBROUTINE CAL_EDGELENGTH
      SUBROUTINE CAL_JAC( Jac, LON, LAT )
        IMPLICIT NONE
        REAL*8:: Jac, LON(3), LAT(3)
        REAL*8:: XR(2), XS(2)
        XR = 0.5d0*(/ LON(2) - LON(1), LAT(2) - LAT(1) /) ; 
        XS = 0.5d0*(/ LON(3) - LON(1), LAT(3) - LAT(1) /) ; 
        Jac = (XR(1)*XS(2) - XR(2)*XS(1)) ; ! Find Jacobian
        RETURN ;
      END SUBROUTINE CAL_JAC
!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into the cylindrical
!    *
!    mapping coordinates.
!    *
!    Lon,Lat must be in radians.
!    *
!
!    CPB: I pulled this from ADCIRC and modified it to always use
!    Mercator for this purpose
!                                                                             *
!******************************************************************************
      SUBROUTINE CYLINDERMAP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0
     
      X= Rearth*(RLAMBDA-RLAMBDA0)*COS(PHI0)  ;
      Y= Rearth*log(tan(PHI) + 1.0d0/cos(PHI))*cos(PHI0) ;          
      
      RETURN ;
      END SUBROUTINE CYLINDERMAP
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE COMPUTE_CYLINPROJ_SFAC( SLAMV, SFEAV )
!----------------------------------------------------------------------
!----------------------------------------------------------------------
        IMPLICIT NONE
! dummy !
        REAL*8,INTENT(IN),DIMENSION(:) :: SLAMV, SFEAV
        REAL*8,ALLOCATABLE,DIMENSION(:) :: SFMX, SFMY
        REAL*8 :: SFCT, SFCX, SFCY, TANPHI, YCSFAC
       
        ! local !
        INTEGER :: I, MEXP
        REAL*8 :: RFAC1, RFAC2, RFAC3

!.... Set:
!.... SFCT, SFCX, SFCY, SFMX, SFMY, YCSFAC
      ALLOCATE( SFMX(NP), SFMY(NP) )
!.... accordingly in order to determine the forms of equations to be
!solved
!     
!.....DEFAULT:

      MEXP = MOD(22,20) ; ! = 2 if mercator
!......ICS = 20 : Equal area projection
!......    = 21 : CPP
!.....     = 22 : Mercartor
!......from WP
      DO I = 1, NP
         SFCT = cos(SFEAV(I))**MEXP ;           
         SFCX = cos(SFEA0)*( cos(SFEAV(I))**(MEXP - 1) ) ;
         SFCY = cos(SFEA0)**(MEXP - 1) ; 
         
         SFMX(I) = cos(SFEA0)/COS(SFEAV(I)) ;
         SFMY(I) = SFMX(I)**(MEXP - 1) ;
         TANPHI = (TAN(SFEAV(I))/Rearth) ;
         YCSFAC = cos( SFEAV(I) ) ; 
      END DO
        
!     COMPUTE ELEMENT AVERAGE FROM NODAL VECTORS ADJUSTING EQUATIONS
!     TO CYLINDER COORDINATES
      CALL SFAC_ELEAVG( SFMXEle, SFMX, NM, NE ) ; 
      CALL SFAC_ELEAVG( SFMYEle, SFMY, NM, NE ) ;         
      DEALLOCATE( SFMY, SFMX )
!      Only need the  averages of these... 
        
      RETURN ; 
!----------------------------------------------------------------------
      END SUBROUTINE COMPUTE_CYLINPROJ_SFAC
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE SFAC_ELEAVG( SFAELES, SFNODES, NM, NE ) 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NM(:,:), NE
      REAL*8, dimension(:),INTENT(IN) :: SFNODES
      REAL*8,DIMENSION(:),INTENT(INOUT) :: SFAELES
      
      ! local
      INTEGER:: IE, NMI1, NMI2, NMI3
      REAL*8 :: ONETRD 
      
      ONETRD = 1d0/3d0 ;
      DO IE = 1, NE
         NMI1 = NM(IE,1) ; 
         NMI2 = NM(IE,2) ; 
         NMI3 = NM(IE,3) ;  
         SFAELES(IE) = (SFNODES(NMI1) + &
         SFNODES(NMI2) + SFNODES(NMI3))*ONETRD ; 
      END DO
      
      RETURN ;
!----------------------------------------------------------------------
      END SUBROUTINE SFAC_ELEAVG 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE Write_Check_Of_Mesh(filename,outvalue)
      IMPLICIT NONE
      integer :: ii
      CHARACTER(len=12) :: filename
      REAL*8,DIMENSION(:) :: outvalue
      OPEN(15,FILE=trim(filename),STATUS='replace', &
          access='sequential',&
          action='write')
      write(15,*) 'test'
      write(15,*) NE, NP
      do ii = 1,NP
         WRITE(15,1000) ii, slam(ii), sfea(ii), outvalue(ii)
      end do
      do ii = 1,NE
         WRITE(15,*) ii, '3', NM(ii,1), nm(ii,2), nm(ii,3)
      end do
 1000 FORMAT(I8,3(E20.12))
      write(15,*) '0'
      write(15,*) '0'
      write(15,*) '0'
      write(15,*) '0'
      write(15,*) '0'
      CLOSE(15)
      END SUBROUTINE write_check_of_mesh
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE DfDxy_nodal(ConnTable, NumEle, NumNodes, Fv, Dphi1, &
                             Dphi2, Dphi3, EleAreas, Tot_Area, DFDxy)
!----------------------------------------------------------------------
!     This subroutine calculates:
!        DF_j/Dx = \SUM_{i=1}^3 ( F_i Dphi_i )
!     Inputs:
!        - ConnTable = connectivity table of mesh
!        - NumEle    = number of elements in mesh
!        - NumNodes  = number of nodes in mesh
!        - Fv        = functional values at nodes in mesh (size 1xNp)
!        - Dphii     = ai/2Aj OR bi/2Aj depending on if this is a
!                    y-derivative or x-derivative. i denotes node
!                    number, j denotes element. 
!        - EleAreas  = areas of all elements (size 1xNumEle)
!        
!     Outputs:
!        - DfDxy     = area averaged derivatives of F_i on a nodal basis.
!                    = (1/TotalArea_i)*\SUM_{j=1}^{NumNei} ( A_j DF_j/dx )
!----------------------------------------------------------------------
      IMPLICIT NONE
      ! inputs
      INTEGER,INTENT(IN),DIMENSION(:,:) :: ConnTable
      INTEGER,INTENT(IN) :: NumEle, NumNodes
      REAL*8,INTENT(IN),DIMENSION(:) :: Fv, Dphi1, Dphi2, Dphi3, &
                                        EleAreas, Tot_Area
      ! outputs
      REAL*8,INTENT(INOUT),DIMENSION(:) :: DFDxy
      ! local
      integer :: ii, NM1, NM2, NM3
      REAL*8  :: DFDxy_ele      
      ! loop through elements and find elemental derivative. Then area
      ! weight and put on the nodes of the element
      write(6,*) 'inside dfdx calculation'
      DO ii = 1,NumEle
         NM1 = ConnTable(ii,1)
         NM2 = ConnTable(ii,2)
         NM3 = ConnTable(ii,3)
         DFDxy_ele = Dphi1(ii)*Fv(NM1) + &
                     Dphi2(ii)*FV(NM2) + &
                     Dphi3(ii)*FV(NM3)
         DFDxy(NM1) = DFDxy(NM1) + DFDxy_ele*EleAreas(ii)
         DFDxy(NM2) = DFDxy(NM2) + DFDxy_ele*EleAreas(ii)
         DFDxy(NM3) = DFDxy(NM3) + DFDxy_ele*EleAreas(ii)
      ENDDO
      ! divide by total area
      DO ii = 1,NumNodes
         DFDxy(ii) = DFDxy(ii)/Tot_Area(ii)
      ENDDO

      END SUBROUTINE DfDxy_nodal
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE Get_Interpolation_Weights()
!----------------------------------------------------------------------
!     This subroutine finds the interpolation weights and indices to 
!     linearly interpolate between the GOFS3.1 grid and ADCIRC grid. It
!     uses the mesh info (SLAM, SFEA, NP, DP) and populates the indices
!     (indxy, indz) based on the GOFS information (BC3D_Lon, BC3D_Lat,
!     and BC3D_Z). To determine the xy information it uses bl_interp
!     from ADCIRC.
!----------------------------------------------------------------------
      IMPLICIT NONE
      ! logical for first call
      LOGICAL,SAVE :: FirstCall=.TRUE.
      ! looping variable
      INTEGER :: ii, iz
      ! lon, lat, and depth of node
      REAL*8 :: xx, yy, bb
      ! indices and weights for each loop iteration
      INTEGER :: indt(4)
      REAL*8 :: wt(4)
      ! allocate on first call
      write(6,*) "getting interp weights"
      IF (FirstCall) THEN
         FirstCall = .FALSE.
         ALLOCATE( INDXY(4,NP), WEIGHTS(4,NP), INDZ(NP) )
      ENDIF
      ! lets get those indices
      DO ii = 1,NP
         xx = SLAM(ii)
         IF (MINVAL(BC3D_Lon).GE.0d0.AND.xx.LT.0d0) xx = xx + 360d0
         yy = SFEA(ii)
         bb = DP(ii)
         ! xy first
         CALL bl_interp(NX,BC3D_Lon,NY,BC3D_Lat,xx,yy,indt,wt)
         INDXY(:,ii)   = indt
         WEIGHTS(:,ii) = wt
         ! z next
         DO iz = 1,NZ
            IF (BC3D_Z(iz).LT.bb) INDZ(ii) = iz
         ENDDO
      ENDDO
      END SUBROUTINE Get_Interpolation_Weights
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE bl_interp(xp, x_array, yp, y_array, x, y, ii, w)
      ! This function uses bilinear interpolation to get the
      ! interpolation weights, w and the indices for interpolation
      ! Assumed to be sampled on a regular grid, with the grid x
      ! values
      ! specified by x_array and the grid y values specified by
      ! y_array
      implicit none
      integer, intent(in) :: xp, yp           
      real*8, dimension(xp), intent(in) :: x_array
      real*8, dimension(yp), intent(in) :: y_array
      real*8, intent(in) :: x,y
      real*8, intent(out) :: w(4)
      integer, intent(out) :: ii(4)
      real*8 :: denom, x1, x2, y1, y2
      real*8 :: x2x1, y2y1, x2x, y2y, xx1, yy1
      integer  :: i, j, ir, jr
      i = binarysearch(xp, x_array, x)
      j = binarysearch(yp, y_array, y)
      ! Make sure we start from one
      if (i.eq.0) then
         i = 1; ir = 1
      else
         ir = i + 1
         if (ir > xp) then
            if (3d0*x_array(1) - 2d0*x_array(2) + &
               360d0 < x_array(xp)) then
               ! Wrap the longitude around
               ir = 1
            else
               ir = xp
            endif
         endif
      endif
      x1 = x_array(i)
      x2 = x_array(ir)
          
      if (j.eq.0) then
         j = 1; jr = 1
      else
         jr = j + 1
         if (jr > yp) jr = yp
      endif
      y1 = y_array(j)
      y2 = y_array(jr)
      ! Lon, lat to meters 
      if (ir.eq.i) then
         x2x1 = 1.0d0
         x2x  = 0.0d0
         xx1  = 1.0d0
      else
         x2x1 = haversine(x1,x2,y,y)
         x2x  = haversine(x,x2,y,y)
         xx1  = haversine(x1,x,y,y)
      endif
      if (jr.eq.j) then
         y2y1 = 1.0d0
         y2y  = 0.0d0
         yy1  = 1.0d0
      else
         y2y1 = haversine(x,x,y1,y2)
         y2y  = haversine(x,x,y,y2)
         yy1  = haversine(x,x,y1,y)
      endif 
      denom = x2x1*y2y1
         
      if (denom < 1e-12) then
         write(6,*) 'Denominator is small',denom
         write(6,*) x,y,x1,y1,x2,y2,i,j,ir,jr
      endif
        
      ! Output the indices
      ii   = [i, j, ir, jr]
      ! Output the weights
      w(1) = x2x*y2y/denom
      w(2) = xx1*y2y/denom
      w(3) = x2x*yy1/denom
      w(4) = xx1*yy1/denom
      END SUBROUTINE bl_interp
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
      FUNCTION binarysearch(length, array, value, delta) 
          ! Given an array and a value, returns the index of the element
          ! that
          ! is closest to, but less than, the given value.
          ! Uses a binary search algorithm.
          ! "delta" is the tolerance used to determine if two values are
          ! equal
          ! if ( abs(x1 - x2) <= delta) then
          !    assume x1 = x2
          ! endif
          implicit none
          integer, intent(in) :: length
          real*8, dimension(length), intent(in) :: array
          !f2py depend(length) array
          real*8, intent(in) :: value
          real*8, intent(in), optional :: delta
          integer :: binarysearch
          integer :: left, middle, right, orientation
          real*8 :: d
          if (present(delta) .eqv. .true.) then
             d = delta
          else
             d = 1d-9
          endif
          orientation = 1
          if (array(2) < array(1)) orientation = -1
       
          left = 1
          right = length
          do
             if (left > right) exit
             middle = nint((left+right) / 2.0d0)
             if ( abs(array(middle) - value) <= d) then
                binarySearch = middle
                return
             endif
             select case(orientation)
             case(1) 
                if (array(middle) > value) then
                   right = middle - 1
                else
                   left = middle + 1
                end if
             case(-1)
                if (array(middle) < value) then
                   right = middle - 1
                else
                   left = middle + 1
                end if
             end select
          end do
          binarysearch = right
      END FUNCTION binarysearch
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE GETOGCMFILE(flag,TimeIn,FileOut)
!----------------------------------------------------------------------
!     This subroutine returns the name of the file we want for a
!     particular input time (timein) and flag. If flag = 1 then we
!     want the temperature and salinity file. If flag = 2 then we want
!     the velocity file. It returns a string "fileout"
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: FLAG
      TYPE(DATETIME),INTENT(IN) :: TimeIn
      CHARACTER(LEN=80),INTENT(OUT) :: FileOut
      INTEGER :: kk
      ! check the bounds first
      IF ( TIMEIN.LT.FN%TSNAP(1) .OR.&
           TIMEIN.GT.FN%TSNAP(FN%NSNAPS) ) THEN
         WRITE(6,*) 'ERROR: This time is not available.'
         CALL ABORT
      ENDIF
      DO kk = 1,FN%NSNAPS-1
         IF ( TIMEIN.GE.FN%TSNAP(kk) .AND.&
              TIMEIN.LT.FN%TSNAP(kk+1) ) THEN
            EXIT
         ENDIF
      ENDDO
      IF ( FLAG.EQ.1 ) THEN
         FileOut = FN%TSFILE(kk)
      ELSEIF ( FLAG.EQ.2 ) THEN
         FileOut = FN%UVFILE(kk)
      ENDIF
!----------------------------------------------------------------------
      END SUBROUTINE GETOGCMFILE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE READOGCMDATA(OGCMFILE)
!----------------------------------------------------------------------
!     This subroutine reads the text file "OGCMFILE" that contains 
!     the timesnaps and filenames that we want. It is structured as:
!
!             nTimeSnaps
!             TSNAP1, TSFILE1, UVFILE1, SSHFILE1
!             TSNAP2, TSFILE2, UVFILE2, SSHFILE2
!             ...
!             TSNAPnTimeSnaps, TSFILEnTimeSnaps,...
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=16),INTENT(IN) :: OGCMFILE
      INTEGER :: nsnaps, kk
      CHARACTER(len=16) :: tsnap ! dummy string b4 datetime conv
      CHARACTER(LEN=40) :: file1,file2,file3
      OPEN(12,file=OGCMFILE,status='old')
      READ(12,*) nsnaps
      FN%NSNAPS=nsnaps
      ALLOCATE( FN%TSNAP(nsnaps), FN%TSFILE(nsnaps),&
                FN%UVFILE(nsnaps), FN%SSHFILE(nsnaps) )
      DO kk = 1,nsnaps
         READ(12,*) tsnap, file1, file2, file3
         FN%TSNAP(kk) = STRPTIME(trim(tsnap),"%Y-%m-%d %H:%M")
         FN%TSFILE(kk) = trim(file1)
         FN%UVFILE(kk) = trim(file2)
         FN%SSHFILE(kk) = trim(file3)
      ENDDO
      CLOSE(12)
!----------------------------------------------------------------------
      END SUBROUTINE READOGCMDATA
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END PROGRAM OGCM_DL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
