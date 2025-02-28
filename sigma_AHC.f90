  subroutine sigma_SOHC
  !This is a second order hall conductance.
  !\begin{equation}
  !\begin{aligned}
  !   &\chi_{x y y}^{\mathrm{int}}=\int_{\mathrm{BZ}} \frac{d \boldsymbol{k}}{(2 \pi)^2} \Lambda_{x y y}(\boldsymbol{k}),\\
		&\text { with }\\
		&\Lambda_{x y y}(\boldsymbol{k})=-\sum_n\left(v_x G_{y y}^n-v_y G_{x y}^n\right) \partial_{E_n} f_0,
  !\end{aligned}
  !\end{equation}

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: m, n, ie
     integer :: ierr, knv3

     real(dp) :: kdotr, mu, Beta_fake
     real(dp) :: k(3)

     real(dp) :: time_start, time_end,epsilon

     ! eigen value of H
	  real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :)
     complex(dp), allocatable :: UU_dag(:, :)

     !> velocities
     complex(dp), allocatable :: vx(:, :), vy(:, :), vz(:, :)
     complex(dp), allocatable :: vxmat(:, :), vymat(:, :), vzmat(:, :)
    
     !> Berry curvature
     complex(dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
     complex(dp), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi
     real(dp)  :: diffFermi

     allocate( W (Num_wann))
     allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann), vz(Num_wann, Num_wann))
     allocate( vxmat(Num_wann, Num_wann), vymat(Num_wann, Num_wann), vzmat(Num_wann, Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( Amat(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( UU_dag(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum))
     allocate( sigma_tensor_ahc_mpi(3, OmegaNum))
     allocate(Omega_x(Num_wann), Omega_y(Num_wann), Omega_z(Num_wann))
     allocate(Omega_x_t(Num_wann), Omega_y_t(Num_wann), Omega_z_t(Num_wann))
     sigma_tensor_ahc    = 0d0
     sigma_tensor_ahc_mpi= 0d0
     vx=0d0
     vy=0d0
     vz=0d0
     Hamk_bulk=0d0
     Amat= 0d0
     UU_dag=0d0
     UU= 0d0
     
     !> energy
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie

     ! Check OmegaNum and OmegaMin/OmegaMax 
     if (OmegaNum < 1 .or. OmegaMin > OmegaMax) then
         write(*,*) "Error: OmegaNum or Omega range invalid."
         stop
     endif

     knv3= Nk1*Nk2*Nk3

     call now(time_start) 
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
           call now(time_end) 
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_latticegauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
       !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)
        if (any(W /= W)) then
            write(*,*) "Error: Eigenvalues contain NaN."
            stop
        endif
  
        vx= 0d0; vy= 0d0; vz= 0d0; vxmat=0d0; vymat=0d0 ; vzmat=0d0
        UU_dag= conjg(transpose(UU))
        do iR= 1, Nrpts
           if (ndegen(iR) == 0) then
               write(*,*) "Error: ndegen(iR) = 0."
               stop
           endif
           kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
           vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        enddo ! iR
   
        !> unitility rotate velocity
        call mat_mul(Num_wann, vx, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vxmat) 
        call mat_mul(Num_wann, vy, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vymat) 
        call mat_mul(Num_wann, vz, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vzmat) 
   
        Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
        epsilon = 1e-8

        do m= 1, Num_wann
           do n= 1, Num_wann
              if (m==n) cycle
              Omega_x(m)= Omega_x(m)+ vxmat(n, m)*vymat(m,n)/((W(m)-W(n))**3+epsilon**3)
              Omega_y(m)= Omega_y(m)+ vymat(n, m)*vymat(m,n)/((W(m)-W(n))**3+epsilon**3)
           enddo ! m           
        enddo ! n
   
        Omega_x= 2.d0*Real(Omega_x)
        Omega_y= 2.d0*Real(Omega_y)

        if (any(Omega_x /= Omega_x) .or. any(Omega_y /= Omega_y)) then
            write(*,*) "Error: Berry curvature contains NaN."
            stop
        endif  
 
        !> consider the Fermi-distribution according to the brodening Earc_eta
        Beta_fake= 1d0/Eta_Arc

        do ie=1, OmegaNum
           mu = energy(ie)
!          if (abs(Beta_fake * (W(m) - mu)) > 1d2) then
!             diffFermi = 0d0
!          else
!             diffFermi = -Beta_fake / (Exp(Beta_fake * (W(m) - mu)) + 1d0) /(Exp(-Beta_fake * (W(m) - mu)) + 1d0)
!          endif
           do m= 1, Num_wann
              diffFermi=-Beta_fake/(Exp(Beta_fake*(W(m)-mu))+1d0)/(Exp(-Beta_fake*(W(m)-mu))+1d0)

              Omega_x_t(m)= Omega_x(m)*diffFermi*vy(m,m)
              Omega_y_t(m)= Omega_y(m)*diffFermi*vx(m,m)
           enddo
           sigma_tensor_ahc_mpi(1, ie)= sigma_tensor_ahc_mpi(1, ie)+ sum(Omega_x_t)
           sigma_tensor_ahc_mpi(2, ie)= sigma_tensor_ahc_mpi(2, ie)+ sum(Omega_y_t)
           sigma_tensor_ahc_mpi(3, ie)= sigma_tensor_ahc_mpi(3, ie)+ sum(Omega_x_t-Omega_y_t)
        enddo ! ie
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif
     sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume   !/CellVolume*24341d0  ! in (Omega*cm)^-1

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(unit=outfileindex, file='sigma_ahc.txt')
        write(outfileindex, '("#",a)')' Anomalous hall conductivity in unit of (Omega*cm)^-1'
        write(outfileindex, '("#",a13, 20a16)')'Eenergy (eV)', '\sigma_xy', '\sigma_yz', '\sigma_zx'
        do ie=1, OmegaNum
           write(outfileindex, '(200E16.8)')energy(ie), sigma_tensor_ahc(3, ie), &
                                                        sigma_tensor_ahc(1, ie), &
                                                        sigma_tensor_ahc(2, ie)

        enddo
        close(outfileindex)
     endif

     deallocate( W , vx, vy, vz, Hamk_bulk, Amat, UU, UU_dag, energy)
     deallocate( sigma_tensor_ahc, sigma_tensor_ahc_mpi)
 
     return
  end subroutine sigma_SOHC
