  subroutine Berry_curvature_singlek_numoccupied_slab_total(k, Omega_z)

     use wmpi
     use para
     implicit none

     !> input parameter, k is in unit of reciprocal lattice vectors
     real(dp), intent(in) :: k(2)
     complex(dp), intent(out) :: Omega_z(1)

     integer :: m, n, Mdim, i1, i2
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_slab(:, :)
     complex(dp), allocatable :: vx(:, :), vy(:, :)
     complex(dp), allocatable :: Vij_x(:, :, :), Vij_y(:, :, :)

     !> leading order of the hamiltonian matrix for slab system specified with Nslab
     Mdim = Num_wann*Nslab

     allocate(W(Mdim))
     allocate(vx(Mdim, Mdim), vy(Mdim, Mdim))
     allocate(UU(Mdim, Mdim), UU_dag(Mdim, Mdim), Hamk_slab(Mdim, Mdim))
     allocate(Amat(Mdim, Mdim), DHDk(Mdim, Mdim, 3), DHDkdag(Mdim, Mdim, 3))
     allocate(Vij_x(-ijmax:ijmax, Num_wann, Num_wann))
     allocate(Vij_y(-ijmax:ijmax, Num_wann, Num_wann))
     W=0d0; vx= 0d0; vy= 0d0; UU= 0d0; UU_dag= 0d0

     call ham_qlayer2qlayer_velocity(k, Vij_x, Vij_y) 
     do i1=1, nslab
        ! i2 row index
        do i2=1, nslab
          if (abs(i2-i1).le.ijmax)then
            vx((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Vij_x(i2-i1,1:Num_wann,1:Num_wann)
            vy((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Vij_y(i2-i1,1:Num_wann,1:Num_wann)
          endif 
        enddo ! i2
     enddo ! i1

     ! calculation slab hamiltonian by a direct Fourier transformation of HmnR
     call ham_slab(k, Hamk_slab)

     !> diagonalization by call zheev in lapack
     UU=Hamk_slab
     call eigensystem_c( 'V', 'U', Mdim, UU, W)
    !call zhpevx_pack(hamk_slab,Mdim, W, UU)

     UU_dag= conjg(transpose(UU))
 
     !> unitility rotate velocity
     call mat_mul(Mdim, vx, UU, Amat) 
     call mat_mul(Mdim, UU_dag, Amat, vx) 
     call mat_mul(Mdim, vy, UU, Amat) 
     call mat_mul(Mdim, UU_dag, Amat, vy) 

     Omega_z=0d0
     do m= 1, NumOccupied*Nslab
        do n= NumOccupied*Nslab+1, Mdim
           Omega_z(1)= Omega_z(1)+ vy(n, m)*vy(m, n)/((W(m)-W(n))**3)
        enddo ! m
!        Omega_z(1)= Omega_z(1)* vx(n,m)
     enddo ! n

     Omega_z= REAL(Omega_z)/Angstrom2atomic**2
!     Omega_z= -aimag(Omega_z*2d0)/Angstrom2atomic**2

!     do m= 1, NumOccupied*Nslab
!        do n= NumOccupied*Nslab+1, Mdim
!           Omega_z= Omega_z*vy(n,m)
!           print *, vx(n,m)
!!           Omega_z(1)= Omega_z(1)+ vx(n,m)*vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
!        enddo ! m
!     enddo ! n
!     Omega_z = Omega_z*sum(vy(1:NumOccupied*Nslab, 1:NumOccupied*Nslab))

     return
  end subroutine Berry_curvature_singlek_numoccupied_slab_total

  subroutine Berry_curvature_slab
     !> Calculate Berry curvature 
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> Aug. 06 2018 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, ierr, Mdim, i, j,i1,i2,n,m

     real(dp) :: k(2)  

     !> k points slice
     real(dp), allocatable :: k12(:, :)
     real(dp), allocatable :: k12_xyz(:, :)
   
     real(dp), external :: norm

     !> Berry curvature  (k)
     complex(dp), allocatable :: Omega_z(:)
     complex(dp), allocatable :: Omega(:)
     complex(dp), allocatable :: Omega_mpi(:)
     complex(dp), allocatable :: vx(:, :), vy(:, :)
     complex(dp), allocatable :: Vij_x(:, :, :), Vij_y(:, :, :)


     Mdim = Num_wann*Nslab

     allocate( k12(2, Nk1*Nk2))
     allocate( k12_xyz(2, Nk1*Nk2))
     allocate( Omega_z(Mdim))
     allocate( Omega    (Nk1*Nk2))
     allocate( Omega_mpi(Nk1*Nk2))
     allocate(vx(Mdim, Mdim), vy(Mdim, Mdim))
     allocate(Vij_x(-ijmax:ijmax, Num_wann, Num_wann))
     allocate(Vij_y(-ijmax:ijmax, Num_wann, Num_wann))

     k12=0d0
     k12_xyz=0d0
     omega= 0d0
     omega_mpi= 0d0
     vx= 0d0; vy= 0d0    
     !> k12 is centered at K2d_start
     ik=0
     do i= 1, nk1
        do j= 1, nk2
           ik=ik+1
           k12(:, ik)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)- (K2D_vec1+K2D_vec2)/2d0
           k12_xyz(:, ik)= k12(1, ik)* Ka2+ k12(2, ik)* Kb2
        enddo
     enddo

     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0) write(stdout, *)'Berry curvature ik, nk1*nk2 ', ik, Nk1*Nk2

        !> diagonalize hamiltonian
        k= k12(:, ik)

        Omega_z= 0d0

        call Berry_curvature_singlek_numoccupied_slab_total(k, Omega_z(1))

        Omega(ik) = sum(Omega_z)

     enddo  !ik

     Omega_mpi= 0d0

#if defined (MPI)
     call mpi_allreduce(Omega,Omega_mpi,size(Omega_mpi),&
                       mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
     Omega_mpi= Omega
#endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_slab.dat')
        write(outfileindex, '(20a28)')'# Unit of Berry curvature is Angstrom^2'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', &
           'Omega_z'
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')k12_xyz(:, ik)*Angstrom2atomic, real(Omega_mpi(ik))/Angstrom2atomic**2
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
    if (cpuid==0) then

        open(unit=outfileindex, file='Berrycurvature_slab.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"set output 'QM.eps'"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"set border lw 1"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set xtics"
        write(outfileindex, '(a)')"set ytics"
        write(outfileindex, '(a)')"set xlabel 'K_1 (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'K_2 (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel offset 1, 0"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"set palette defined (0 'blue', 0.5 'white')"
        write(outfileindex, '(a)')'#set palette defined (-10 '#194eff', 0 'white', 10 'red' )'
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"splot 'Berrycurvature_slab.dat' u 1:2:3 w pm3d"
        close(outfileindex)

endif



#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( k12)
     deallocate( k12_xyz)
     deallocate( Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine Berry_curvature_slab
