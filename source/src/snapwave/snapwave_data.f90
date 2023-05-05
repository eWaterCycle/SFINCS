   module snapwave_data
   !
   real*8,  dimension(:),       allocatable    :: x, y                    ! x,y coordinates of grid (double precision)
   real*4,  dimension(:),       allocatable    :: xs, ys                  ! x,y coordinates of grid (single precision)
   real*8                                      :: xmn, ymn
   integer*1, dimension(:),     allocatable    :: msk                     ! mask array (0=land, 1=inner point, 2=boundary)
   real*4,  dimension(:),       allocatable    :: zb, depth               ! bed level, water depth
   real*4,  dimension(:),       allocatable    :: dhdx, dhdy              ! water depth gradients
   real*4,  dimension(:),       allocatable    :: kwav, nwav              ! wave number, ratio Cg/C
   real*4,  dimension(:),       allocatable    :: kwav_ig, nwav_ig        ! wave number, ratio Cg/C
   real*4,  dimension(:),       allocatable    :: C, Cg                   ! wave celerity, group velocity0
   real*4,  dimension(:),       allocatable    :: C_ig, Cg_ig             ! wave celerity, group velocity0
   real*4,  dimension(:),       allocatable    :: fw                      ! friction coefficient. Double for interpolation routine
   real*4,  dimension(:),       allocatable    :: fw_ig                   ! friction coefficient
   real*4,  dimension(:),       allocatable    :: H, H_ig                 ! rms wave height
   real*4,  dimension(:),       allocatable    :: Tp                      ! peak wave period
   real*4,  dimension(:),       allocatable    :: Dw,Df                   ! dissipation due to breaking, bed friction
   real*4,  dimension(:),       allocatable    :: Sw,St                   ! Source term for windwave growth, wave energy, wave period
   real*4,  dimension(:,:),     allocatable    :: see                     ! source term windwave growth, directionally distributed
   real*4,  dimension(:),       allocatable    :: F                       ! wave force Dw/C/rho/depth
   real*4,  dimension(:),       allocatable    :: Fx, Fy                  ! wave force Dw/C/rho/depth
   real*4,  dimension(:),       allocatable    :: u10, u10dir             ! wind speed and wind direction
   real*4                                      :: u10dmean                ! average wind direction
   real*4,  dimension(:),       allocatable    :: thetam                  ! mean wave direction
   real*4                                      :: thetamean               ! mean wave direction
   real*4,  dimension(:),       allocatable    :: buf                     ! buffer for writing output to netcdf
   real*4,  dimension(:),       allocatable    :: buf1                    ! buffer for writing output to netcdf
   real*4,  dimension(:,:),     allocatable    :: buf2                    ! buffer for writing output to netcdf
   integer, dimension(:,:),     allocatable    :: kp                      ! surrounding points for each grid point (unstructured)
   logical, dimension(:),       allocatable    :: inner                   ! mask for inner points (not on any boundary)
   integer, dimension(:),       allocatable    :: neumannconnected        ! neumann boundary indices connected to inner points
   integer                                     :: noneumannpts            ! number of neumann boundary points
   real*4,  dimension(:),       allocatable    :: theta                   ! wave angles,sine and cosine of wave angles
   real*4,  dimension(:),       allocatable    :: dist                    ! relative distribution of energy over theta bins
   real*4,  dimension(:),       allocatable    :: theta360                ! wave angles,sine and cosine of wave angles
   integer, dimension(:),       allocatable    :: i360                    ! reference between partial thea grid and full 360 deg theta grid
   real*4,  dimension(:,:),     allocatable    :: windspread360           ! wind input distribution array full 360 deg theta grid
   real*4,  dimension(:,:),     allocatable    :: windspread              ! wind input distribution array
   integer, dimension(:,:,:),   allocatable    :: prev                    ! two upwind grid points per grid point and wave direction
   integer, dimension(:,:,:),   allocatable    :: prev360                 ! two upwind grid points per grid point and wave direction
   integer, dimension(:,:),   allocatable      :: prevu10                 ! two upwind grid points for wind dir per grid point 
   real*4,  dimension(:,:,:),   allocatable    :: w                       ! weights of upwind grid points, 2 per grid point and per wave direction
   real*4,  dimension(:,:,:),   allocatable    :: w360                    ! weights of upwind grid points, 2 per grid point and per wave direction
   real*4,  dimension(:,:),   allocatable      :: wu10                    ! weights of upwind grid points for wind dir, 2 per grid point    
   real*4,  dimension(:,:),     allocatable    :: ds                      ! distance to interpolated upwind point, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ds360                   ! distance to interpolated upwind point, per grid point and direction
   real*4,  dimension(:),       allocatable    :: dsu10                   ! distance to interpolated upwind point for wind dir, per grid point 
   real*4,  dimension(:,:),     allocatable    :: ctheta                  ! refraction speed, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ctheta_ig               ! refraction speed, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ctheta360               ! refraction speed, per grid point and direction
   !   real*4,  dimension(:),       allocatable    :: xn,yn,zn               ! coordinates of nodes of unstructured grid
   real*4,  dimension(:),       allocatable    :: dzdx,dzdy               ! bed slopes at nodes of unstructured grid
   integer, dimension(:,:),     allocatable    :: face_nodes             ! node numbers connected to each cell
   integer, dimension(:,:),     allocatable    :: edge_nodes             ! node numbers connected to each edge
   real*4,  dimension(:),       allocatable    :: bndindx
   real*4,  dimension(:),       allocatable    :: tau
   !   real*4,  dimension(:,:),     allocatable    :: Fluxtab
   !   real*4,  dimension(:),       allocatable    :: F0tab
   real*4                                      :: Hmax
   real*4,  dimension(:),       allocatable    :: sinhkh
   real*4,  dimension(:),       allocatable    :: Hmx
   real*4,  dimension(:),       allocatable    :: sinhkh_ig
   real*4,  dimension(:),       allocatable    :: Hmx_ig
   real*4,  dimension(:,:),     allocatable    :: ee                      ! directional energy density
   real*4,  dimension(:,:),     allocatable    :: ee_ig                   ! directional infragravity energy density
   !
   integer*4,  dimension(:),     allocatable    :: index_snapwave_in_quadtree
   integer*4,  dimension(:),     allocatable    :: index_quadtree_in_snapwave
   
   character*256 :: trefstr_iso8601
   character*41  :: treftimefews
   character*15  :: trefstr
   character*15  :: tstartstr
   character*15  :: tstopstr
   !
   ! Boundary conditions (single point, time-varying)
   !
   character*232                               :: jonswapfile             ! filename of time-varying wave and water level data
   character*232                               :: windlistfile            ! filename of time-varying wind magnitude and direction
   
   !
   ! Boundary conditions (space- and time-varying)
   !
   character*232                               :: bndfile
   character*232                               :: encfile
   character*232                               :: bhsfile
   character*232                               :: btpfile
   character*232                               :: bwdfile
   character*232                               :: bdsfile
   character*232                               :: bzsfile

   integer                                     :: nwbnd                   ! number of support points wave boundary
   integer                                     :: ntwbnd                  ! number of time points wave boundary
   real*4                                      :: tpmean_bwv              ! mean tp over boundary points for given time
   real*4                                      :: wdmean_bwv              ! mean wave direction for given time, used to make theta grid
   real*4                                      :: zsmean_bwv              ! mean water level for given time, used to make theta grid
   real*8,  dimension(:),     allocatable      :: x_bwv                   ! x coordinates of boundary points
   real*8,  dimension(:),     allocatable      :: y_bwv                   ! y coordinates of boundary points
   real*4,  dimension(:),     allocatable      :: t_bwv                   ! times (s) of wave boundary conditions
   integer                                     :: n_bndenc                ! number of boundary enclosure points
   real*8,  dimension(:),     allocatable      :: x_bndenc                ! x coordinates of boundary enclosure points
   real*8,  dimension(:),     allocatable      :: y_bndenc                ! y coordinates of boundary enclosure points
   integer                                     :: n_neu                   ! number of neumann polyline points
   real*8,  dimension(:),     allocatable      :: x_neu                   ! x coordinates of neumann polyline
   real*8,  dimension(:),     allocatable      :: y_neu                   ! y coordinates of neumann polyline
   integer                                     :: nnmb                    ! number of neumann grid points
   real*4,  dimension(:),     allocatable      :: hst_bwv                 ! wave height at boundary points for given time
   real*4,  dimension(:),     allocatable      :: tpt_bwv                 ! wave period at boundary points for given time
   real*4,  dimension(:),     allocatable      :: wdt_bwv                 ! wave direction at boundary points for given time
   real*4,  dimension(:),     allocatable      :: dst_bwv                 ! directional spreading at boundary points for given time
   real*4,  dimension(:),     allocatable      :: zst_bwv                 ! water level at boundary points for given time
   real*4,  dimension(:,:),     allocatable    :: eet_bwv                 ! directional spectra at boundary points for given time
   !
   real*4,  dimension(:,:),     allocatable    :: hs_bwv                  ! wave height for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: tp_bwv                  ! wave period for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: wd_bwv                  ! wave direction (nautical deg) for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: ds_bwv                  ! directional spreading (deg) for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: zs_bwv                  ! water level for all boundary locations and time points
   !
   integer                                     :: ntu10bnd                ! number of time points wind forcing boundary conditions
   integer*4                                   :: itwindbndlast
   real*4,  dimension(:),     allocatable      :: t_u10_bwv               ! times (s) of wind boundary conditions
   real*4,  dimension(:,:),     allocatable    :: u10_bwv                 ! wind speed for all grid cells and time points
   real*4,  dimension(:,:),     allocatable    :: u10dir_bwv              ! wind direction for all grid cells and time points   
   !real*4,  dimension(:),     allocatable      :: u10t_bwv                ! wind speed at all grid cells for given time
   !real*4,  dimension(:),     allocatable      :: u10dirt_bwv             ! wave direction at all grid cells for given time
     
   !
   integer*4,          dimension(:),   allocatable :: index_v_m
   integer*4,          dimension(:),   allocatable :: index_v_n
   integer*4,          dimension(:),   allocatable :: index_v_nm
   integer*4,          dimension(:,:), allocatable :: index_g_nm
   !
   integer*4                                       :: itwbndlast
   integer*4,          dimension(:),   allocatable :: nmindbnd            ! index of grid point at grid boundary
   integer*4,          dimension(:),   allocatable :: neubnd              ! index of grid point at neumann grid boundary
   integer*4,          dimension(:),   allocatable :: ind1_bwv_cst        ! index to closest wave boundary point
   integer*4,          dimension(:),   allocatable :: ind2_bwv_cst        ! index to second closest wave boundary point
   real*4,             dimension(:),   allocatable :: fac_bwv_cst         ! weight of closest point

   integer*4 mmax
   integer*4 nmax
   real*4 dx
   real*4 dy
   real*4 x0
   real*4 y0
   real*4 rotation
   real*4 cosrot
   real*4 sinrot
   !
   ! Timing
   !
   real*8 :: tstart
   real*8 :: tstop
   real*4 :: timestep
   !
   character*232                               :: map_filename
   character*232                               :: his_filename
   !
   ! Local input variables
   !
   !   integer                                   :: nx,ny           ! number of grid cells in two directions
   !   integer                                   :: m,n             ! number of grid points in two directions
   real*4                                    :: dtheta              ! theta grid resolution
   real*4                                    :: sector              ! theta grid sector
   real*4                                    :: fw0                 ! uniform wave friction factor
   real*4                                    :: fw0_ig              ! uniform wave friction factor (ig waves)
   real*4                                    :: Tpini               ! initial condition for the wave period
   real*4                                    :: zsini               ! initial condition for the water level
   logical                                   :: broadcast_Tp_bwc    ! switch that initializes domain with average boundary wave period
   real*4                                    :: fwcutoff            ! depth below which to apply space-varying fw
   real*4                                    :: snapwave_alpha,gamma         ! coefficients in Baldock breaking dissipation model
   real*4                                    :: hmin                ! minimum water depth
   character*232                             :: gridfile            ! name of gridfile (Delft3D .grd format)
   integer                                   :: sferic              ! sferical (1) or cartesian (0) grid
   integer                                   :: niter               ! maximum number of iterations
   real*4                                    :: crit                ! relative accuracy for stopping criterion
   character*232                             :: depfile             ! name of bathymetry file (Delft3D .dep format)
   character*232                             :: fwfile              ! name of bed friction factor file (Delft3D .dep format)
   character*232                             :: upwfile             ! name of upwind neighbors file
   character*232                             :: mskfile             ! name of mask file
   character*232                             :: indfile             ! name of index file
   character*232                             :: obsfile             ! name of observation points file
   character*256                             :: fwstr, fw_igstr     ! name of bed friction files OR uniform value of fw and fw_ig
   character*256                             :: u10str, u10dirstr   ! name of bed friction files OR uniform value of u10 and u10dir   
   integer                                   :: nobs
   real*8,    dimension(:),    allocatable   :: xobs
   real*8,    dimension(:),    allocatable   :: yobs
   integer ,  dimension(:,:),  allocatable   :: irefobs
   integer ,  dimension(:),    allocatable   :: nrefobs
   real*8  ,  dimension(:,:),  allocatable   :: wobs
   character*32, dimension(:), allocatable   :: nameobs         ! names of observation points
   real*4,    dimension(:),    allocatable   :: hm0obs
   real*4,    dimension(:),    allocatable   :: hm0igobs
   real*4,    dimension(:),    allocatable   :: hm0xobs
   real*4,    dimension(:),    allocatable   :: hm0yobs
   real*4,    dimension(:),    allocatable   :: tpobs
   real*4,    dimension(:),    allocatable   :: wdobs
   real*4                                    :: dt              ! time step (no limitation)
   real*4                                    :: tol             ! tolerance(m) for boundary points
   integer                                   :: no_nodes        ! number of unstructured grid nodes
   integer                                   :: no_faces        ! number of unstructured grid cells
   integer                                   :: no_edges        ! number of unstructured grid edges
   !   integer                                   :: k,j,k1,k2
   !   integer                                   :: irec,nrec
   integer                                   :: ntab = 10
   integer                                   :: nHrel
   !   real*4                                    :: xmin,xmax       ! minimum and maximum x coordinates
   !   real*4                                    :: ymin,ymax       ! minimum and maximum y coordinates
   !   real*4                                    :: ddir,alfa
   character*232                             :: Htabname,Dwtabname,Ftabname,Cgtabname,cthetafactabname,hhtabname
   character*232                             :: neumannfile
   real*4, dimension(:), allocatable         :: xb,yb,xneu,yneu

   !
   character*3                               :: outputformat
   !
   ! Local constants
   !
   real*4                                    :: rho             ! water density
   real*4                                    :: pi              ! water density
   real*4                                    :: g               ! acceleration of gravity
   real*4                                    :: t0,t1,t2,t3,t4  ! timers
   integer                                   :: nb
   integer                                   :: np
   integer                                   :: ntheta
   integer                                   :: ntheta360
   !
   logical                                   :: igwaves
   !real*4                                    :: crit
   integer                                   :: nr_sweeps
   !
   logical                                   :: restart
   logical                                   :: coupled_to_sfincs
   !
   integer                                   :: nr_quadtree_points
   !
   ! Output variables
   !
   integer                                   :: map_Hm0
   integer                                   :: map_Hig
   integer                                   :: map_Tp
   integer                                   :: map_dir
   integer                                   :: map_cg
   integer                                   :: map_Dw
   integer                                   :: map_Df
   integer                                   :: map_Sw
   integer                                   :: map_St
   integer                                   :: map_ee
   !
   end module snapwave_data
