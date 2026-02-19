#include "tlab_error.h"

module SpecialForcing
    use TLab_Constants, only: wp, wi, pi_wp, tag_scal
    use TLab_Constants, only: efile, lfile, MAX_VARS, MAX_PARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Grid, only: x, y, z
    use Tlab_Background, only: qbg, sbg
    use Profiles, only: Profiles_Calculate
    implicit none
    private

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt
    type(term_dt), public, protected :: forcingProps              ! Forcing parameters

    public :: SpecialForcing_Initialize
    public :: SpecialForcing_Source
    public :: SpecialForcing_Sinusoidal, SpecialForcing_Sinusoidal_NoSlip   ! Taylor-Green vortex
    public :: GravityWave_Polarization

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_HOMOGENEOUS = 1
    integer, parameter :: TYPE_SINUSOIDAL = 2
    integer, parameter :: TYPE_SINUSOIDAL_NOSLIP = 3
    integer, parameter :: TYPE_RAND_MULTIPLICATIVE = 4
    integer, parameter :: TYPE_RAND_ADDIVTIVE = 5
    integer, parameter :: TYPE_WAVEMAKER = 6
    integer, parameter :: TYPE_WAVEMAKERTANH = 7

    integer, parameter :: nwaves_max = 3                ! maximum number of waves in wavemaker
    integer :: nwaves                                   ! number of waves in wavemaker
    real(wp) amplitude(4, nwaves_max)                   ! wave amplitudes in x, y, z + buoyancy
    real(wp) wavenumber(3, nwaves_max)                  ! wavenumbers in x, y, z
    real(wp) frequency(nwaves_max)                      ! wave frequencies
    real(wp) envelope(6)                                ! (x,y,z) position, size

    real(wp), allocatable, target :: tmp_envelope(:, :, :) ! arrays for this routine
    real(wp), allocatable, target :: tmp_phase(:, :, :)

contains
    !########################################################################
    !########################################################################
    subroutine SpecialForcing_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax, inb_scal
        use IO_Fields, only: IO_Read_Fields
        use TLab_Time, only: itime
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_offset_i, ims_offset_j
#endif
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy, i, j, k, iwave, idsp, jdsp
        real(wp) :: dummy(MAX_PARS), rx, ry, rz
        real(wp) :: params(0)

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'SpecialForcing'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<values>')
        call TLab_Write_ASCII(bakfile, '#Wave#=<amplitude,wavenumber,angle,frequency>')
        call TLab_Write_ASCII(bakfile, '#Envelope=<x,y,z,size>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; forcingProps%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'homogeneous') then; forcingProps%type = TYPE_HOMOGENEOUS; 
        elseif (trim(adjustl(sRes)) == 'random') then; forcingProps%type = TYPE_RAND_MULTIPLICATIVE; 
        elseif (trim(adjustl(sRes)) == 'sinusoidal') then; forcingProps%type = TYPE_SINUSOIDAL; 
        elseif (trim(adjustl(sRes)) == 'wavemaker') then; forcingProps%type = TYPE_WAVEMAKER; 
        elseif (trim(adjustl(sRes)) == 'wavemakertanh') then; forcingProps%type = TYPE_WAVEMAKERTANH; 
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in SpecialForcing.Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        forcingProps%vector(:) = 0.0_wp; forcingProps%active(:) = .false.
        if (forcingProps%type /= TYPE_NONE) then
            call ScanFile_Char(bakfile, inifile, block, 'Vector', '1.0, 0.0, 0.0', sRes)
            idummy = 3
            call LIST_REAL(sRes, idummy, forcingProps%vector)

            if (abs(forcingProps%vector(1)) > 0.0_wp) then; forcingProps%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Ox.'); end if
            if (abs(forcingProps%vector(2)) > 0.0_wp) then; forcingProps%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Oy.'); end if
            if (abs(forcingProps%vector(3)) > 0.0_wp) then; forcingProps%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Forcing along Oz.'); end if

            forcingProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '1.0, 1.0, 0.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, forcingProps%parameters)

            select case (forcingProps%type)
            case (TYPE_HOMOGENEOUS)

            case (TYPE_WAVEMAKER, TYPE_WAVEMAKERTANH)
                forcingProps%active(1:3) = .true.       ! default is active in x, y, z momentum equations

                do nwaves = 1, nwaves_max
                    write (sRes, *) nwaves
                    call ScanFile_Char(bakfile, inifile, block, 'Wave'//trim(adjustl(sRes)), 'void', sRes)
                    if (trim(adjustl(sRes)) /= 'void') then
                        idummy = 4
                        call LIST_REAL(sRes, idummy, dummy)
                        if (idummy /= 4) then
                            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in '//trim(adjustl(block))//'.Wave.')
                            call TLab_Stop(DNS_ERROR_OPTION)
                        end if
                        dummy(3) = dummy(3)*pi_wp/180._wp                   ! from degree to radians
                        wavenumber(1, nwaves) = dummy(2)*cos(dummy(3))      ! x-wavenumber
                        wavenumber(2, nwaves) = 0.0_wp
                        wavenumber(3, nwaves) = dummy(2)*sin(dummy(3))      ! z-wavenumber
                        amplitude(1, nwaves) = dummy(1)*sin(dummy(3))       ! in x equation
                        amplitude(2, nwaves) = 0.0_wp
                        amplitude(3, nwaves) = -dummy(1)*cos(dummy(3))      ! in z equation
                        frequency(nwaves) = dummy(4)
                        amplitude(4, nwaves) = (-sbg(1)%delta/sbg(1)%thick)**2/frequency(nwaves)*amplitude(3,nwaves) ! sbg not yet initialized
                    else
                        exit
                    end if
                end do
                nwaves = nwaves - 1                                         ! correct for the increment in the loop

                ! forcingProps%active(2) = .false.                            ! only active in x and z
            end select

        end if

        !########################################################################
        !               Initialize the envelope for the wave maker
        !########################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; jdsp = ims_offset_j
#else
        idsp = 0; jdsp = 0
#endif

        select case (forcingProps%type)
        case (TYPE_WAVEMAKER, TYPE_WAVEMAKERTANH)
            envelope(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Envelope', '1.0, 1.0, 1.0, 1.0, 1.0, 1.0', sRes) ! position and size
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, envelope)

            call TLab_Allocate_Real(__FILE__, tmp_envelope, [imax, jmax, kmax], 'tmp-wave-envelope')
            call TLab_Allocate_Real(__FILE__, tmp_phase, [imax, kmax, nwaves], 'tmp-wave-phase')

            ! Add a random perturbation to the envelope in order to support turbulence onset
            ! forcingProps%parameters(2) = perturbation amplitude
            if (forcingProps%parameters(2) < 1.0_wp .AND. forcingProps%parameters(2) > 0.0_wp) then
                call IO_Read_Fields(trim(adjustl(tag_scal))//'rand', imax, jmax, kmax, itime, inb_scal, 1, tmp_envelope, params)
                ! Normalize random field
                tmp_envelope = tmp_envelope/max(abs(maxval(tmp_envelope)), abs(minval(tmp_envelope)))
                tmp_envelope = tmp_envelope*forcingProps%parameters(2)
            end if
        end select

        select case (forcingProps%type)
        case (TYPE_WAVEMAKER)
            ! envelope(1:3) are postions in x, y, and z

            envelope(4) = abs(envelope(4)) ! make sure the size parameter is positive
            envelope(5) = abs(envelope(5)) ! make sure the size parameter is positive
            envelope(6) = abs(envelope(6)) ! make sure the size parameter is positive

            dummy(1) = 0.5_wp/envelope(4)**2.0_wp ! width in x-direction
            dummy(2) = 0.5_wp/envelope(5)**2.0_wp ! width y-direction
            dummy(3) = 0.5_wp/envelope(6)**2.0_wp ! width z-direction

            do k = 1, kmax
                do j = 1, jmax
                    do i = 1, imax
                        rx = x%nodes(idsp + i) - envelope(1)
                        ry = y%nodes(jdsp + j) - envelope(2)
                        rz = z%nodes(k) - envelope(3)
                        tmp_envelope(i, j, k) = (tmp_envelope(i, j, k) + 1.0_wp)*exp(- dummy(1)*rx*rx - dummy(2)*ry*ry - dummy(3)*rz*rz)
                        do iwave = 1, nwaves
                            tmp_phase(i, k, iwave) = rx*wavenumber(1, iwave) + rz*wavenumber(3, iwave)
                        end do
                    end do
                end do
            end do

            ! ! Get buoyancy freq. profile and store into wrk1d; Needed for b-amplitude
            ! call fdm_der1_Z%compute(nx*ny, Profiles_Calculate(sbg(is), z%nodes(k)), wrk1d)

        case (TYPE_WAVEMAKERTANH)
            ! Tanh profile for wave maker envelope only a function of z
            ! envelope has only two parameters: inflection point position and width
            ! Tanh profile thus only for simulating a horizontally plane wave
            
            ! envelope(1) is position of inflection point
            ! envelope(2) is profile width

            do k = 1, kmax
                do j = 1, jmax
                    do i = 1, imax
                        rx = x%nodes(idsp + i) - envelope(1)
                        ry = y%nodes(jdsp + j) - envelope(2)
                        rz = z%nodes(k) - envelope(3)
                        tmp_envelope(i, j, k) = (tmp_envelope(i, j, k) + 1.0_wp)*0.5*(1.0 + tanh((z%nodes(k) - envelope(1))/(2.0*envelope(2))))
                        do iwave = 1, nwaves
                            tmp_phase(i, k, iwave) = rx*wavenumber(1, iwave) + rz*wavenumber(3, iwave)
                        end do
                    end do
                end do
            end do
            
        end select

        return
    end subroutine SpecialForcing_Initialize

!########################################################################
!########################################################################
    subroutine SpecialForcing_Source(locProps, nx, ny, nz, iq, time, q, h, tmp)
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, iq
        real(wp), intent(in) :: time
        real(wp), intent(in) :: q(nx, ny, nz)
        real(wp), intent(inout) :: h(nx, ny, nz)
        real(wp), intent(inout) :: tmp(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) iwave, i, k
        ! real(wp) shift

        !########################################################################
        select case (locProps%type)

        case (TYPE_HOMOGENEOUS)
            tmp = locProps%parameters(1)

        case (TYPE_RAND_MULTIPLICATIVE)
            call random_number(tmp)

            tmp = (tmp*2.0_wp - 1.0_wp)*locProps%parameters(1)
            tmp = tmp*h

        case (TYPE_SINUSOIDAL)
            tmp = locProps%parameters(1)*sin(locProps%parameters(2)*time)

        case (TYPE_SINUSOIDAL_NOSLIP)

        case (TYPE_WAVEMAKER, TYPE_WAVEMAKERTANH)
            do k = 1, nz
                do i = 1, nx
                    do iwave = 1, nwaves
                        tmp(i, 1:ny, k) = amplitude(iq, iwave)
                        tmp(i, 1:ny, k) = tmp(i, 1:ny, k)*sin(tmp_phase(i, k, iwave) - frequency(iwave)*time)
                    end do
                end do
                tmp(1:nx, 1:ny, k) = Profiles_Calculate(qbg(iq), z%nodes(k)) + tmp(1:nx, 1:ny, k)
            end do
            tmp = (tmp - q)*tmp_envelope*locProps%parameters(1)
            tmp = tmp*ft(time)
        end select

        return
    end subroutine SpecialForcing_Source


    ! Perhaps one could avoid this routine by giving tmp_phase a shift in source above?
    ! tmp_phase(b) = tmp_phase(w) + pi/2
    subroutine GravityWave_Polarization(locProps, nx, ny, nz, is, time, s, tmp)
        ! Buoyancy part of a GW is phase shifted by pi/2 compared to the 
        ! velocity components. We impose a wave with u=Usin(phase), w=Wsin(phase), 
        ! where U=Asin(Theta) and W=-Acos(Theta) (incompressibility).
        ! Linearized Boussinesq => Fourier space => B=-iN^2/freq*W
        ! => b=-N^2/freq*W*cos(phase)
        ! Amplitude is set by assuming a linear buoyancy profile (N(z)=N=const).

        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        real(wp), intent(in) :: time
        real(wp), intent(in) :: s(nx, ny, nz)
        real(wp), intent(inout) :: tmp(nx, ny, nz)
        integer(wi) iwave, i, k

        select case (forcingProps%type)
        case (TYPE_WAVEMAKER, TYPE_WAVEMAKERTANH)
            do k = 1, nz
                do i = 1, nx
                    do iwave = 1, nwaves
                        tmp(i, 1:ny, k) = amplitude(3+is, iwave)
                        tmp(i, 1:ny, k) = -tmp(i, 1:ny, k)*cos(tmp_phase(i, k, iwave) - frequency(iwave)*time)      ! Wave signal
                    end do
                end do
                tmp(1:nx, 1:ny, k) = Profiles_Calculate(sbg(is), z%nodes(k)) + tmp(1:nx, 1:ny, k)                 ! Add Background
            end do
            tmp = (tmp - s)*tmp_envelope*locProps%parameters(1)
            tmp = tmp*ft(time)
        end select
 
        return
    end subroutine
    
    
    real(wp) function ft(t) 
        real(wp), intent(in) :: t
        ft = tanh(0.1_wp*t)
    end function ft

    !########################################################################
    ! Sinusoidal forcing; Taylor-Green vortex
    !########################################################################
    subroutine SpecialForcing_Sinusoidal(nx, ny, nz, time, visc, u, v, h_u, h_v)
        integer(wi) nx, ny, nz
        real(wp) time, visc
        real(wp), dimension(nx, ny, nz) :: u, v, h_u, h_v

        ! -----------------------------------------------------------------------
        real(wp) omega, sigma, amplitude
        integer(wi) i, k

        !########################################################################
        omega = 2.0_wp*pi_wp
        sigma = 2.0_wp*omega*omega*visc

        amplitude = -(1.0_wp + (sigma/omega)**2)*omega
        amplitude = amplitude*sin(omega*time)
        do k = 1, nz
            do i = 1, nx
                u(i, :, k) = sin(x%nodes(i)*omega)*cos(z%nodes(k)*omega)
                v(i, :, k) = -cos(x%nodes(i)*omega)*sin(z%nodes(k)*omega)
            end do
        end do

        amplitude = sin(omega*time)

        h_u = h_u + amplitude*u
        h_v = h_v + amplitude*v

        return
    end subroutine SpecialForcing_Sinusoidal

    !########################################################################
    ! Velocity field with no-slip
    !########################################################################
    subroutine SpecialForcing_Sinusoidal_NoSlip(nx, ny, nz, time, visc, h1, h2, tmp1, tmp2, tmp3, tmp4)
        use OPR_Partial, only: OPR_Partial_X, OPR_Partial_Z, OPR_P1, OPR_P2_P1

        integer(wi) nx, ny, nz
        real(wp) time, visc
        real(wp), dimension(nx, ny, nz) :: h1, h2
        real(wp), dimension(nx, ny, nz) :: tmp1, tmp2, tmp3, tmp4

        ! -----------------------------------------------------------------------
        integer(wi) i, k

        ! #######################################################################
        do k = 1, nz
            do i = 1, nx
                !     tmp1(i, :, k) = sin(2.0_wp*pi_wp*x%nodes(i))*       sin(C_4_R*pi_wp*z%nodes(j))
                !     tmp2(i, :, k) =-cos(2.0_wp*pi_wp*x%nodes(i))*(1.0_wp-cos(C_4_R*pi_wp*z%nodes(j)))*0.5_wp
                tmp1(i, :, k) = sin(pi_wp*x%nodes(i))*sin(pi_wp*x%nodes(i))*sin(2.0_wp*pi_wp*z%nodes(k))
                tmp2(i, :, k) = -sin(2.0_wp*pi_wp*x%nodes(i))*sin(pi_wp*z%nodes(k))*sin(pi_wp*z%nodes(k))
            end do
        end do

        ! Time terms
        h1 = h1 - tmp1*2.0_wp*pi_wp*sin(2.0_wp*pi_wp*time)
        h2 = h2 - tmp2*2.0_wp*pi_wp*sin(2.0_wp*pi_wp*time)

        ! velocities
        tmp1 = tmp1*cos(2.0_wp*pi_wp*time)
        tmp2 = tmp2*cos(2.0_wp*pi_wp*time)

        ! Diffusion and convection terms in Ox momentum eqn
        call OPR_Partial_Z(OPR_P2_P1, nx, ny, nz, tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_Partial_X(OPR_P2_p1, nx, ny, nz, tmp1, tmp4, tmp3)
        h1 = h1 - visc*(tmp4) + (tmp3*tmp1)

        ! Diffusion and convection terms in Oy momentum eqn
        call OPR_Partial_Z(OPR_P2_P1, nx, ny, nz, tmp2, tmp4, tmp3)
        h2 = h2 - visc*(tmp4) + (tmp3*tmp2)

        call OPR_Partial_X(OPR_P2_p1, nx, ny, nz, tmp2, tmp4, tmp3)
        h2 = h2 - visc*(tmp4) + (tmp3*tmp1)

        ! #######################################################################
        do k = 1, nz
            do i = 1, nx
                !     tmp1(i, :, k) = cos(C_4_R*pi_wp*x%nodes(i))*(2.0_wp-cos(C_4_R*pi_wp*z%nodes(j)))/C_8_R &
                !          - 0.5_wp*(sin(2.0_wp*pi_wp*z%nodes(j)))**4
                tmp1(i, :, k) = sin(2.0_wp*pi_wp*x%nodes(i))*sin(2.0_wp*pi_wp*z%nodes(k))
                tmp1(i, :, k) = tmp1(i, :, k)*(cos(2.0_wp*pi_wp*time))**2
            end do
        end do

        ! Pressure gradient
        call OPR_Partial_X(OPR_P1, nx, ny, nz, tmp1, tmp2)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, tmp1, tmp3)

        h1 = h1 + tmp2
        h2 = h2 + tmp3

        return
    end subroutine SpecialForcing_Sinusoidal_NoSlip

end module SpecialForcing
