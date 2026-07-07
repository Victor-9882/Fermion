PROGRAM SPLINE_2D_FERMION_BOSON
  IMPLICIT DOUBLE PRECISION (a-h, o-z)

  PARAMETER (NMG=6, NMZ=6, NMA=2*NMG*NMZ, LWORK=10*NMA)

  DOUBLE PRECISION :: XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA)
  DOUBLE PRECISION :: VR(NMA,NMA), VL(1,1)
  DOUBLE PRECISION :: ALPHAR(NMA), ALPHAI(NMA), BETA(NMA)
  DOUBLE PRECISION :: WR(NMA), WI(NMA), WORK(LWORK)
  DOUBLE PRECISION :: zv(NMZ+1), gv(NMG+1)
  DOUBLE PRECISION :: c(2, NMG,NMZ), XG(NMA), YG(NMA)
  DOUBLE PRECISION :: Xq(1000), DXq(1000), Yp(1000), DYp(1000)
  DOUBLE PRECISION :: Wv(1000), DWv(1000)
  DOUBLE PRECISION :: Uu(1000), DUu(1000)

  INTEGER :: ii, i, j, NPARAM, INFO, iw, f, s, t, r !f vai numerar a equação acoplada
  INTEGER :: N_intervalG, N_intervalZ, NCOL
  INTEGER :: Nnz(100), Nng(100), Nnv(100), Nnu(100)
  INTEGER :: Nz, Ng, Nv, Nu
  DOUBLE PRECISION :: e, mf, ms, Mtot, mu, kappa, PI, gam0, Alfa, m, xi, u
  DOUBLE PRECISION, EXTERNAL :: KERNEL_LOWER, KERNEL_UPPER
  CHARACTER :: kernel_type

  open(unit=10, file="autovalores.dat",  STATUS="UNKNOWN")
  open(unit=11, file="autovetores.dat",  STATUS="UNKNOWN")
  open(unit=12, file="alfa.dat",         STATUS="UNKNOWN")
  open(unit=13, file="ZMatrix",          STATUS="UNKNOWN")
  open(unit=14, file="XMatrix",          STATUS="UNKNOWN")
  open(unit=15, file="outro.dat",         STATUS="UNKNOWN")
  open(unit=16, file="coeficientes.dat", STATUS="UNKNOWN")
  open(unit=20, file="inputs.dat",       STATUS="UNKNOWN")

  e  =  0.000001d0
  PI = DACOS(-1.d0)

  READ(20,*) NPARAM
  READ(20,*) kernel_type
  DO i = 1, NPARAM
    READ(20,*) Nnz(i), Nng(i), Nnv(i), Nnu(i)
  END DO
  CLOSE(20)

  Mtot  = 1.90d0
  mf     = 1.d0
  ms = 1.d0
  m = (ms + mf)/2
  mu    = 0.d0
  kappa = SQRT(m**2 - 0.25d0*Mtot**2)
  gam0  = 3.d0
  xi = 0.9999d0

  WRITE(10,'(A,I0,A,I0,A,I0)') "NMA: ",NMA," NMG: ",NMG," NMZ: ",NMZ
  WRITE(10,'(A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10)') &
    "Mtot: ",Mtot," mf: ",mf," ms: ",ms," m: ",m," mu: ",mu," kappa: ",kappa," e: ",e," gam0: ",gam0," xi: ",xi
  WRITE(10,*) "--------------------------------------------------"

  WRITE(12,'(A,I0,A,I0,A,I0)') "NMA: ",NMA," NMG: ",NMG," NMZ: ",NMZ
  WRITE(12,'(A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10)') &
    "Mtot: ",Mtot," mf: ",mf," ms: ",ms," m: ",m," mu: ",mu," kappa: ",kappa," e: ",e," gam0: ",gam0," xi: ",xi
  WRITE(12,*) "--------------------------------------------------"

  iw          = 14
  N_intervalZ = (NMZ-1)/2
  N_intervalG = (NMG-1)/2
  NCOL        = 2

  DO ii = 1, NPARAM
    print*, ii

    Nz = Nnz(ii)
    Ng = Nng(ii)
    Nv = Nnv(ii)
    Nu = Nnu(ii)

    ! z-mesh via collocação
    call G1D(iw, -1.d0, N_intervalZ, 1.d0, 1.d0, Xq)
    call COLLOC(iw, NCOL, N_intervalZ, Xq, XG)
    do i = 1, 2*N_intervalZ
      zv(i+1) = XG(i)
    end do
    zv(1)   = -0.999999999d0
    zv(NMZ) =  0.99999999d0

    ! gamma-mesh via collocação
    call G1D(iw, 0.d0, N_intervalG, 1.d0, gam0, Yp)
    call COLLOC(iw, NCOL, N_intervalG, Yp, YG)
    do i = 1, 2*N_intervalG
      gv(i+1) = YG(i)
    end do
    gv(1)   = 0.d0
    gv(NMG) = gam0

    if (xi /= 1.d0) gv(1) = 0.001d0

    call SPLGR1(zv, NMZ)
    call SPLGR2(gv, NMG)

    ! Nós e pesos de Gauss-Legendre
    CALL legauss(0.d0, 1.d0, Nz, Xq,  DXq,  1.d-15)
    CALL legauss(0.d0, gam0, Ng, Yp,  DYp,  1.d-15)
    CALL legauss(0.d0, 1.d0, Nv, Wv,  DWv,  1.d-15)
    CALL legauss(0.d0, 1.d0, Nu, Uu,  DUu,  1.d-15)


    !WRITE(15,'(A,I0,A,I0)') "Nv: ",Nv,"  Nu: ",Nu
    !WRITE(15,'(A20,2X,A20,2X,A20)') "v", "Uu(t)", "u = Uu(t)*(1-v)"
   ! do r = 1, Nv
     ! do t = 1, Nu
        !WRITE(15,'(F20.12,2X,F20.12,2X,F20.12)') Wv(r), Uu(t), Uu(t)*(1.d0 - Wv(r))
      !end do
    !end do
    !WRITE(15,*) ""


    if (xi == 1.d0 .or. kernel_type /= 'v') then
      CALL BUILD_ZMATRIX(zv, gv, NMZ, NMG, NMA, m, mu, kappa, PI, &
                         Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                         ms, mf, Mtot, ZMATRIX, XMATRIX, kernel_type, xi, &
                         Nu, Uu, DUu)
    else
      CALL BUILD_ZMATRIX_CALIBRE(zv, gv, NMZ, NMG, NMA, m, mu, kappa, PI, &
                                  Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                                  ms, mf, Mtot, ZMATRIX, XMATRIX, kernel_type, xi, &
                                  Nu, Uu, DUu)
    end if

    do i = 1, NMA
      XMATRIX(i,i) = XMATRIX(i,i) + e
    end do

    DO i = 1, NMA
      WRITE(13,'(9999ES16.8)') (ZMATRIX(i,j), j=1, NMA)
    END DO

    DO i = 1, NMA
      WRITE(14,'(9999ES16.8)') (XMATRIX(i,j), j=1, NMA)
    END DO

    CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
               ALPHAR, ALPHAI, BETA, VL, 1, VR, NMA, &
               WORK, LWORK, INFO)

    IF (INFO /= 0) THEN
      PRINT *, 'ERRO NO DGGEV: INFO = ', INFO
      STOP
    END IF

    DO i = 1, NMA
      if (ABS(BETA(i)) > 1.d-16) then
        WR(i) = ALPHAR(i) / BETA(i)
        WI(i) = ALPHAI(i) / BETA(i)
      else
        WR(i) =  1.d+16
        WI(i) =  0.d0
      end if
    END DO

    WRITE(10,'(A,I0,A,I0,A,I0,A)') "autovalores_Nz",Nz,"_Ng",Ng,"_Nv",Nv,".dat"
    WRITE(12,'(A,I0,A,I0,A,I0,A,I0,A)') &
      "Numericamente.Collocation.autovalores_Nz",Nz,"_Ng",Ng,"_Nv",Nv,"_Nu",Nu,".dat"

    do i = 1, NMA
      WRITE(10,'(I4,2X,F20.12,2X,F20.12)') i, WR(i), WI(i)
    end do

    WRITE(10,*) ""
    Alfa = 1.0d0 / (WR(1) * 8.0d0 * PI)
    WRITE(10,'(A,F20.10)') "Valor de Alfa: ", Alfa
    WRITE(12,'(A,F20.10)') "Valor de Alfa: ", Alfa

    WRITE(10,'(9999ES16.8)') (VR(j,1), j=1, NMA)

    do s = 1, 2
      do j = 1, NMZ
        do i = 1, NMG
          c(s,i,j) = (VR(i + (j-1)*NMG + (s-1)*NMG*NMZ, 1))
        end do
      end do
    end do

    DO s = 1, 2
      DO i = 1, NMG
        WRITE(16,'(9999ES16.8)') (c(s,i,j), j=1, NMZ)
      END DO
    END DO

  END DO

  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(15)
  CLOSE(13)
  CLOSE(14)
  CLOSE(16)

END PROGRAM SPLINE_2D_FERMION_BOSON


! -----------------------------------------------------------------------
! BUILD_ZMATRIX
! Monta ZMATRIX (lado direito) e XMATRIX (lado esquerdo).
! Para trocar o kernel físico, modifique KERNEL_UPPER e KERNEL_LOWER.
! -----------------------------------------------------------------------
SUBROUTINE BUILD_ZMATRIX(zv, gv, Nmz, Nmg, NMA, m, mu, kappa, PI, &
                          Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                          ms, mf, Mtot, ZMATRIX, XMATRIX, kernel_type, xi, &
                          Nu, Uu, DUu)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)

  INTEGER, INTENT(IN) :: Nmz, Nmg, NMA, Nz, Ng, Nv, Nu
  DOUBLE PRECISION, INTENT(IN)  :: zv(*), gv(*)
  DOUBLE PRECISION, INTENT(IN)  :: m, mu, kappa, PI, ms, mf, Mtot, xi
  DOUBLE PRECISION, INTENT(IN)  :: Xq(Nz), DXq(Nz)
  DOUBLE PRECISION, INTENT(IN)  :: Yp(Ng), DYp(Ng)
  DOUBLE PRECISION, INTENT(IN)  :: Wv(Nv), DWv(Nv)
  DOUBLE PRECISION, INTENT(IN)  :: Uu(Nu), DUu(Nu)
  DOUBLE PRECISION, INTENT(OUT) :: ZMATRIX(NMA,NMA), XMATRIX(NMA,NMA)
  CHARACTER, INTENT(IN) :: kernel_type

  DOUBLE PRECISION :: splz_z(Nmz), splz_q(Nmz), splg(Nmg)
  DOUBLE PRECISION :: z, g, gp, zq, v, dzq, dgp, dv
  DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER, KERNEL_LOWER
  INTEGER :: i, j, k, l, p, q, r, index1, index2, f, s

  ZMATRIX = 0.d0
  XMATRIX = 0.d0
  u = 0.d0 !Não depende de u

do s = 1, 2
  do i = 1, Nmg
    g = gv(i)
    print*, i
    do j = 1, Nmz
      z      = zv(j)
      
      call SPLMD1(zv, Nmz, z, splz_z)

      do f = 1, 2
      do k = 1, Nmg
        do l = 1, Nmz
          index1 = (f-1)*Nmg*Nmz + (j-1)*Nmg + i
          index2 = (s-1)*Nmg*Nmz + (l-1)*Nmg + k

          do p = 1, Ng
            gp  = Yp(p)
            dgp = DYp(p)
            call SPLMD2(gv, Nmg, gp, splg)

            if (f == s) then
              ! Lado esquerdo: integra apenas em gp
              XMATRIX(index1,index2) = XMATRIX(index1,index2) + &
                splg(k)*splz_z(l)*dgp / &
                (g + gp + (1.d0 - z**2)*kappa**2 + m**2*z**2)**2
            end if

            do q = 1, Nz

              ! Região superior: z' em [z, 1]
              dzq = (1.d0 - z)*DXq(q)
              zq  = (1.d0 - z)*Xq(q) + z
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)
                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_UPPER(z, zq, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv
              end do

              ! Região inferior: z' em [-1, z]
              dzq = (z + 1.d0)*DXq(q)
              zq  = (z + 1.d0)*Xq(q) - 1.d0
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)
                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_LOWER(z, zq, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv
              end do

            end do  ! q
          end do  ! p
        end do  ! l
      end do  ! k
    end do !f
    end do  ! j
  end do  ! i
end do !s

END SUBROUTINE BUILD_ZMATRIX


SUBROUTINE BUILD_ZMATRIX_CALIBRE(zv, gv, Nmz, Nmg, NMA, m, mu, kappa, PI, &
                          Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                          ms, mf, Mtot, ZMATRIX, XMATRIX, kernel_type, xi, Nu, Uu, DUu)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)

  INTEGER, INTENT(IN) :: Nmz, Nmg, NMA, Nz, Ng, Nv, Nu
  DOUBLE PRECISION, INTENT(IN)  :: zv(*), gv(*)
  DOUBLE PRECISION, INTENT(IN)  :: m, mu, kappa, PI, ms, mf, Mtot, xi
  DOUBLE PRECISION, INTENT(IN)  :: Xq(Nz), DXq(Nz)
  DOUBLE PRECISION, INTENT(IN)  :: Yp(Ng), DYp(Ng)
  DOUBLE PRECISION, INTENT(IN)  :: Wv(Nv), DWv(Nv), Uu(Nu), DUu(Nu)
  DOUBLE PRECISION, INTENT(OUT) :: ZMATRIX(NMA,NMA), XMATRIX(NMA,NMA)
  CHARACTER, INTENT(IN) :: kernel_type

  DOUBLE PRECISION :: splz_z(Nmz), splz_q(Nmz), splg(Nmg)
  DOUBLE PRECISION :: z, g, gp, zq, v, dzq, dgp, dv
  DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER, KERNEL_LOWER
  INTEGER :: i, j, k, l, p, q, r, index1, index2, f, s, t

  ZMATRIX = 0.d0
  XMATRIX = 0.d0

do s = 1, 2
  do i = 1, Nmg
    g = gv(i)
    print*, i
    do j = 1, Nmz
      z      = zv(j)
      
      call SPLMD1(zv, Nmz, z, splz_z)

      do f = 1, 2
      do k = 1, Nmg
        do l = 1, Nmz
         index1 = (f-1)*Nmg*Nmz + (j-1)*Nmg + i
          index2 = (s-1)*Nmg*Nmz + (l-1)*Nmg + k

          do p = 1, Ng
            gp  = Yp(p)
            dgp = DYp(p)
            call SPLMD2(gv, Nmg, gp, splg)

            if (f == s) then
              ! Lado esquerdo: integra apenas em gp
              XMATRIX(index1,index2) = XMATRIX(index1,index2) + &
                splg(k)*splz_z(l)*dgp / &
                (g + gp + (1.d0 - z**2)*kappa**2 + m**2*z**2)**2
            end if

            do q = 1, Nz

              ! Região superior: z' em [z, 1]
              dzq = (1.d0 - z)*DXq(q)
              zq  = (1.d0 - z)*Xq(q) + z
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)

              do t = 1, Nu

                u = Uu(t)*(1.d0-v)
                du = (1.d0 - v)*DUu(t)

                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_UPPER(z, zq, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv*du

              end do
            end do

              ! Região inferior: z' em [-1, z]
              dzq = (z + 1.d0)*DXq(q)
              zq  = (z + 1.d0)*Xq(q) - 1.d0
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)
              do t = 1, Nu
                u = Uu(t)*(1.d0-v)
                du = (1.d0 - v)*DUu(t)
                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_LOWER(z, zq, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv*du
              end do
            end do

            end do  ! q
          end do  ! p
        end do  ! l
      end do  ! k
    end do !f
    end do  ! j
  end do  ! i
end do !s

END SUBROUTINE BUILD_ZMATRIX_CALIBRE


! -----------------------------------------------------------------------
! KERNEL_UPPER  —  integrando para z' em [z, 1]  (eq. 15, theta(z'-z)=1)
! -----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION KERNEL_UPPER(z, zp, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)
  DOUBLE PRECISION, INTENT(IN) :: z, zp, g, gp, v, PI, u
  DOUBLE PRECISION, INTENT(IN) :: m, mu, kappa, ms, mf, Mtot, xi
  INTEGER, INTENT(IN) :: s, f
  CHARACTER, INTENT(IN) :: kernel_type
  DOUBLE PRECISION :: M2_4, ku, lD1, lD2, kD  !começam com letra i-n: precisam de declaração explícita
  LOGICAL, PARAMETER :: DEBUG_COEF = .FALSE.  ! mude para .TRUE. para ativar flags de explosão

  M2_4  = 0.25d0*Mtot**2
  Delta = 0.5d0*(ms - mf)

  ! Eq. (16): D_0(gamma, z)
  D0 = g + (1.d0 - z**2)*kappa**2 + (Delta - m*z)**2

              
  Du_e = v*(1.d0-v)*(zp-z)*( g - (1.d0-z**2)*M2_4 + ms**2 )  &
     + (1.d0+z)*( v*(1.d0-v)*( g + z**2*M2_4 )              &
                + v*(gp + kappa**2)                           &
                + v**2*zp**2*M2_4                             &
                + (1.d0-v)*mu**2 )

  ! Eq. (20): k_u^-
  ku = 0.5d0*Mtot - 2.d0*(g + ms**2)/(Mtot*(1.d0+z))

  Bsf  = 0.d0
  Pij2 = 0.d0
  if (kernel_type == 'v') then
    ! Coeficientes vetoriais de P1
    a     = 1.d0 + 2.d0*(mf)/(Mtot)
    c11_0 = 0.5d0*Mtot**2*(a + 0.5d0*v*zp*a + 0.5d0*(2.d0-v)*z) + (2.d0-v)*g
    c11_1 = -0.5d0*Mtot*(2.d0-v)*(1.d0-z)
    c12_0 = 0.5d0*a*( g*(1.d0-v)*(2.d0+v) + 2.d0*v*(gp+kappa**2)  &
          + 2.d0*(1.d0-v)*mu**2 + M2_4*(v*zp-2.d0)*(z-v*(z-zp)) )  &
          + (g + z*M2_4)*(1.d0-v-0.5d0*zp*v)
    c12_1 = 0.5d0*Mtot*( a*(1.d0-v)*(1.d0+z+0.5d0*(z-zp)*v)        &
          - (1.d0-z)*(1.d0-v-0.5d0*zp*v) )
    c21_0 = 0.5d0*Mtot**2*((2.d0-v)*(2.d0-a) + 2.d0 + zp*v)
    c22_0 = g*(1.d0-v)*(2.d0+v) + 2.d0*v*(gp+kappa**2) + 2.d0*(1.d0-v)*mu**2  &
          + 0.5d0*Mtot**2*( (1.d0-v)*(0.5d0*z*zp*v - z + 2.d0-a)               &
                          + 0.5d0*zp*v*(zp*v - 4.d0+a) )
    c22_1 = 0.5d0*Mtot*(1.d0-v)*(v*(z-zp) + 2.d0*(1.d0+z))

    if (f==1 .and. s==1) Bsf = c11_0 + c11_1*ku
    if (f==1 .and. s==2) Bsf = c12_0 + c12_1*ku
    if (f==2 .and. s==1) Bsf = c21_0
    if (f==2 .and. s==2) Bsf = c22_0 + c22_1*ku

    
  if (xi /= 1.d0) then
! Constantes novas vetoriais de P2 (xi != 1) (Calibre arbitrário)

    f2 = - g * v*(1.d0 - v)*(1.d0 + zp) &
        - v*(1.d0 + z)*(gp + kappa**2) &
        + (Mtot**2 * v*(1.d0 + z) / 4.d0) &
          * (-z*(1.d0 + zp)*(1.d0 - v) + zp*(1.d0 - v - v*zp)) &
        - ms**2 * v*(1.d0 - v)*(zp - z) &
        - mu**2 * (1.d0 + z)*(1.d0 - u - v + xi*u)

    Du_c = (1.d0/4.d0)*v &
                * (4.d0*g + 4.d0*gp - (Mtot**2*(z+1.d0)*(zp+1.d0)*((v-1.d0)*z-v*zp+1.d0)) &
                + mf**2*(z+1.d0) + 2.d0*mf*ms*(z+1.d0) &
                + 4.d0*ms**2*v*z - 4.d0*ms**2*v*zp - 3.d0*ms**2*z + 4.d0*ms**2*zp + ms**2 &
                - 4.d0*g*v - 4.d0*g*v*zp + 4.d0*gp*z + 4.d0*g*zp) &
                - mu**2*(z+1.d0)*(-xi*u+u+v-1.d0)

! Coeficientes vetoriais de P2 (xi != 1) (Calibre arbitrário)

d11_0 = -(v**2 / 32.d0)*(1+z)                                                  &
        * ( -4.d0*ms**2*(-1.d0 + z)                                             &
            + Mtot**2*(1.d0 + z)*(-1.d0 + 2.d0*z - zp)                         &
            - 2.d0*mf*Mtot*(1.d0 + z)*zp                                        &
            + 8.d0*g )                                                           &
        * ( Mtot**2*(1.d0 + z)*(1.d0 + zp)*(-1.d0 + (-1.d0 + v)*z - v*zp)     &
            - 4.d0*ms**2*(-1.d0 + (-2.d0 + v)*z + zp - v*zp)                   &
            + 4.d0*(-1.d0 + v)*(1.d0 + zp)*g )

    d12_0 = v**2*(1+z)*(Mtot**2*(z+1.d0)*(zp+1.d0)*((v-1.d0)*z-v*zp-1.d0) &
     - 4.d0*ms**2*((v-2.d0)*z-v*zp+zp-1.d0) &
     + 4.d0*g*(v-1.d0)*(zp+1.d0)) &
     * (Mtot**3*(z+1.d0)*(z*(v*zp+v+zp-1.d0)-v*zp*(zp+1.d0)) &
     + 2.d0*Mtot**2*mf*(z+1.d0)*((v-1.d0)*z*(zp+1.d0)-zp*(v*zp+v-1.d0)) &
     - 4.d0*Mtot*ms**2*(z*(v+zp-1.d0)-v*zp) &
     + 4.d0*g*Mtot*(v*zp+v+zp-1.d0) &
     - 8.d0*mf*(v-1.d0)*(ms**2*(z-zp)-g*(zp+1.d0)))/(64.d0*Mtot)


    d21_0 = -((1.d0 + z)**2 / 16.d0)                                                  &
        * ( (Mtot*v**2*(2.d0*mf + Mtot*(-1.d0 + zp)))                          &
            * (-Mtot**2*(1.d0 + z)*(1.d0 + zp)*(-1.d0 + (-1.d0 + v)*z - v*zp) &
               + 4.d0*ms**2*(-1.d0 + (-2.d0 + v)*z + zp - v*zp)                &
               - 4.d0*(-1.d0 + v)*(1.d0 + zp)*g) )

d22_0 = (1.d0/32.d0)*v**2*(1+z) &
     * (-(Mtot**2*(z+1.d0)*(zp+1.d0)*((v-1.d0)*z-v*zp-1.d0)) &
     + 4.d0*ms**2*((v-2.d0)*z-v*zp+zp-1.d0) &
     - 4.d0*g*(v-1.d0)*(zp+1.d0)) &
     * (-(Mtot**2*(z+1.d0)*((v-1.d0)*z*(zp+1.d0)-zp*(v*zp+v-2.d0))) &
     + 2.d0*Mtot*mf*(z+1.d0)*zp &
     + 4.d0*(v-1.d0)*(ms**2*(z-zp)-g*(zp+1.d0)))

    d11_2 = v*(1.d0 + z)                                                            &
        * ( -4.d0*ms**2*(-1.d0 + z)                                             &
            + Mtot**2*(1.d0 + z)*(-1.d0 + 2.d0*z - zp)                         &
            - 2.d0*Mtot*mf*(1.d0 + z)*zp                                        &
            + 8.d0*g )                                                           &
        / 8.d0

    d12_2 = (v*(1.d0 + z) / 16.d0)                                                 &
        * ( 4.d0*ms**2*(-1.d0 + zp - 2.d0*v*zp + z*(-3.d0 + 2.d0*v + zp))      &
            - 2.d0*mf*Mtot*(1.d0 + z)*(-1.d0 - 2.d0*(1.d0 + zp)               &
              * (z - v*z + v*zp))                                                &
            - Mtot**2*(1.d0 + z)*(2.d0*z*(-1.d0 + v + v*zp)                    &
              - (1.d0 + zp)*(1.d0 + 2.d0*v*zp))                                 &
            - 8.d0*(-1.d0 + v + v*zp)*g                                           &
            + (8.d0*mf/Mtot)*(ms**2*(-1.d0 - 3.d0*z + 2.d0*v*z                 &
              + 2.d0*zp - 2.d0*v*zp) - 2.d0*(-1.d0 + v)*(1.d0 + zp)*g) )


    d21_2 = -(Mtot*v*(1.d0 + z)**2 / 4.d0)                                         &
        * (2.d0*mf + Mtot*(-1.d0 + zp))


      d22_2 = -(v*(1.d0 + z) / 8.d0)                                                 &
        * ( -2.d0*mf*Mtot*(1.d0 + z)*zp                                         &
            + 4.d0*ms**2*(1.d0 + (3.d0 - 2.d0*v)*z + 2.d0*(-1.d0 + v)*zp)      &
            + Mtot**2*(1.d0 + z)*(-1.d0 + zp - 2.d0*(1.d0 + zp)                 &
              *(z - v*z + v*zp))                                                  &
            + 8.d0*(-1.d0 + v)*(1.d0 + zp)*g )

denom3 = Mtot**4*(1.d0 + z)**4                                             &
       - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)                                 &
       + 16.d0*(ms**2 + g)**2

d11_3 = (Mtot*(z+1.d0)**4 &
     * (Mtot**3*z*(-2.d0*v*(z+zp+2.d0)+3.d0*z+4.d0) &
     - 2.d0*(Mtot/2.d0-(2.d0*(g+ms**2))/(Mtot*(z+1.d0))) &
       * (Mtot**2*(2.d0*v*(z**2+2.d0*z*(zp+2.d0)-zp-2.d0)-3.d0*z**2-8.d0*z+4.d0) &
          + 4.d0*Mtot*mf*(v*(2.d0*z+zp+2.d0)-3.d0*z-2.d0) &
          + 4.d0*g*(3.d0-2.d0*v)) &
     + 4.d0*Mtot**2*mf*z*(v*(zp+2.d0)-2.d0) &
     + 4.d0*Mtot*(2.d0*v-3.d0)*(z-1.d0)*(Mtot/2.d0-(2.d0*(g+ms**2))/(Mtot*(z+1.d0)))**2 &
     - 4.d0*g*Mtot*(2.d0*v*(z+2.d0*zp+6.d0)-3.d0*z-14.d0) &
     + 16.d0*g*mf*(3.d0-2.d0*v))) &
     / (2.d0*(Mtot**4*(z+1.d0)**4-8.d0*Mtot**2*(z+1.d0)**2*(ms**2-g)+16.d0*(g+ms**2)**2))

d12_3 = ((1.d0 + z)**2 / (4.d0*denom3))                                         &
        * ( 2.d0*Mtot**3*mf*(1.d0 + z)**2*(1.d0 + (5.d0 - 2.d0*z)*z             &
            + 4.d0*v**2*(1.d0 + z)**2*zp                                        &
            - 2.d0*v*(1.d0 + 3.d0*z + (3.d0 + z*(4.d0 + 3.d0*z))*zp))           &
            + Mtot**4*(1.d0 + z)**2*(1.d0 + z*(-1.d0 + 10.d0*z)                 &
            + 4.d0*v**2*(1.d0 + z)**2*zp                                        &
            - 2.d0*v*(1.d0 + z + 2.d0*zp + 7.d0*z*zp + z**2*(4.d0 + zp)))       &
            + 32.d0*(mf/Mtot)*(ms**2 + g)*(ms**2*(1.d0 + 3.d0*z + 4.d0*v**2*zp  &
            - 2.d0*v*(1.d0 + z + 3.d0*zp)) + 2.d0*(-1.d0 + v*(-3.d0             &
            + 2.d0*v)*zp)*g)                                                    &
            - 4.d0*Mtot**2*(1.d0 + z)*(ms**2*(2.d0 + z*(-4.d0 + 15.d0*z)        &
            + 2.d0*v*(-2.d0 + z**2*(-5.d0 + zp) + 4.d0*(-1.d0 + v)*zp           &
            + 4.d0*(-2.d0 + v)*zp))                                             &
            + g + 2.d0*(v - 10.d0*z + 8.d0*v*z                                  &
            + v*(7.d0 + 2.d0*z - 4.d0*v*(1.d0 + z)*zp))*g)                      &
            - 8.d0*Mtot*mf*(1.d0 + z)*(ms**2*(2.d0 + (8.d0 - 3.d0*z)*z          &
            + 2.d0*v*(-2.d0 + z*(-4.d0 + z - 4.d0*zp) - 6.d0*zp)                &
            + 8.d0*v**2*(1.d0 + z)*zp)                                          &
            + (-5.d0 + 6.d0*v + 4.d0*z + 4.d0*v*(2.d0 + 3.d0*z                  &
            - 2.d0*v*(1.d0 + z)*zp))*g)                                         &
            + 16.d0*(ms**4*(1.d0 - 3.d0*z + 6.d0*z**2 + 4.d0*v**2*zp            &
            - 2.d0*v*(1.d0 + 2.d0*zp + z*(-1.d0 + 2.d0*z + zp)))                &
            + ms**2*(5.d0 - 15.d0*z + 2.d0*v*(-3.d0 + 5.d0*z                    &
            - (3.d0 - 4.d0*v + z)*zp))*g                                        &
            + 2.d0*(5.d0 + 2.d0*v**2*zp - v*(4.d0 + zp))*g**2) )

d21_3 = -(Mtot**2*(z+1.d0)**4 &
     * (Mtot**2*(-2.d0*v*(z*zp+z-2.d0*zp-4.d0)+z-8.d0) &
     + 2.d0*(Mtot/2.d0-(2.d0*(g+ms**2))/(Mtot*(z+1.d0))) &
       * (Mtot*(2.d0*v*(2.d0*z+zp+1.d0)-6.d0*z-1.d0)+2.d0*mf*(2.d0*v-3.d0)) &
     - 2.d0*Mtot*mf*(2.d0*v*(z+2.d0*zp+4.d0)-3.d0*z-8.d0) &
     + 8.d0*g*(2.d0*v-3.d0))) &
     / (Mtot**4*(z+1.d0)**4-8.d0*Mtot**2*(z+1.d0)**2*(ms**2-g)+16.d0*(g+ms**2)**2)


d22_3 = ((z+1.d0)**2 &
     * (Mtot**4*(z+1.d0)**2*(4.d0*v**2*(z+1.d0)**2*(zp+1.d0) &
        - 2.d0*v*(z**2*(3.d0*zp+4.d0)+z*(3.d0*zp+7.d0)+4.d0*zp+7.d0) &
        + 3.d0*z**2+9.d0*z+10.d0) &
     - 4.d0*Mtot**3*mf*(z+1.d0)**2*(v*((z-1.d0)*zp-2.d0)+z+2.d0) &
     - 4.d0*Mtot**2*(z+1.d0) &
       * (ms**2*(8.d0*v**2*(z+1.d0)*(zp+1.d0)+2.d0*v*(z**2-2.d0*z*(2.d0*zp+5.d0)-7.d0*zp-12.d0)-3.d0*z**2+16.d0*z+16.d0) &
          - g*(8.d0*v**2*(z+1.d0)*(zp+1.d0)-2.d0*v*(6.d0*z*zp+8.d0*z+3.d0*zp+7.d0)+6.d0*z+9.d0)) &
     - 16.d0*Mtot*mf*(z+1.d0)*(g+ms**2*(v*(2.d0*z+zp+2.d0)-3.d0*z-2.d0)+g*v*zp) &
     + 16.d0*(2.d0*v-3.d0)*(g+ms**2)*(ms**2*(2.d0*v*(zp+1.d0)-z-2.d0)+g*(2.d0*v*(zp+1.d0)-1.d0)))) &
     / (2.d0*(Mtot**2*(z+1.d0)**2-4.d0*Mtot*ms*(z+1.d0)+4.d0*(g+ms**2)) &
       * (Mtot**2*(z+1.d0)**2+4.d0*Mtot*ms*(z+1.d0)+4.d0*(g+ms**2)))

d12_4 = (1.d0 + z)**2*(Mtot + 2.d0*mf)                                          &
        / (4.d0*Mtot)

d22_4 = (1.d0 + z)**2 / 2.d0

d11_5 = -(((z+1.0d0)**3*(v*zp+1.0d0)*(Mtot**4*(2.0d0*z**4+3.0d0*z**3+z**2+z+1.0d0) &
        +4.0d0*Mtot**3*mf*z*(z+1.0d0)**2-4.0d0*Mtot**2*(z+1.0d0)*(g+ms**2*(z**2 &
        -2.0d0*z+2.0d0)-4.0d0*g*z)-16.0d0*Mtot*mf*(z+1.0d0)*(ms**2*z-g)-16.0d0 &
        *(g+ms**2)*(ms**2*(z-1.0d0)-2.0d0*g)))/(2.0d0*(Mtot**4*(z+1.0d0)**4 &
        -8.0d0*Mtot**2*(z+1.0d0)**2*(ms**2-g)+16.0d0*(g+ms**2)**2)))

d12_5 = ((1.d0 + z)**3*(1.d0 + v*zp)                                            &
        / (4.d0*(Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2            &
           *(ms**2 - g) + 16.d0*(ms**2 + g)**2)))                                &
        * ( 2.d0*Mtot**3*mf*(1.d0 + z)**2*((-1.d0 + z)*z + v*(1.d0 + z)**2*zp) &
            + Mtot**4*(1.d0 + z)**2*(z - 3.d0*z**2 + v*(1.d0 + z)**2*zp)        &
            - 32.d0*(mf/Mtot)*(ms**2 + g)*(ms**2*(z - v*zp)                     &
              - (1.d0 + v*zp)*g)                                                  &
            + 4.d0*Mtot**2*(1.d0 + z)*(ms**2*(z*(-2.d0 + 5.d0*z)               &
              - 2.d0*v*(1.d0 + z)*zp) + (1.d0 - 6.d0*z                          &
              + 2.d0*v*(1.d0 + z)*zp)*g)                                          &
            + 8.d0*Mtot*mf*(1.d0 + z)*(-ms**2*((-2.d0 + z)*z                    &
              + 2.d0*v*(1.d0 + z)*zp) + (-1.d0 + 2.d0*z                         &
              + 2.d0*v*(1.d0 + z)*zp)*g)                                          &
            + 16.d0*(ms**4*(z - 2.d0*z**2 + v*zp) + ms**2*(-1.d0 + 5.d0*z      &
              + 2.d0*v*zp)*g + (-3.d0 + v*zp)*g**2) )

d21_5 = -(Mtot*(1.d0 + z)**4*(1.d0 + v*zp)                                      &
        / (Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)      &
           + 16.d0*(ms**2 + g)**2))                                               &
        * ( Mtot**3*(1.d0 + z)*(-1.d0 + 3.d0*z)                                 &
            - 2.d0*Mtot**2*mf*(-1.d0 + z**2)                                    &
            - 8.d0*mf*(ms**2 + g)                                                &
            + 4.d0*Mtot*(ms**2*(1.d0 - 2.d0*z) + 3.d0*g) )



     d22_5 = (Mtot*(z+1.d0)**5*(v*zp+1.d0) &
     * (Mtot**3*v*z**2*zp &
     + 2.d0*(Mtot/2.d0-(2.d0*(g+ms**2))/(Mtot*(z+1.d0))) &
       * (-4.d0*g+Mtot**2*z*(2.d0*v*zp+z-2.d0)+4.d0*Mtot*mf*z) &
     - 4.d0*Mtot*(z-v*zp)*(Mtot/2.d0-(2.d0*(g+ms**2))/(Mtot*(z+1.d0)))**2 &
     + 4.d0*g*Mtot*(4.d0*v*zp+z-2.d0)+16.d0*g*mf)) &
     / (2.d0*(Mtot**4*(z+1.d0)**4-8.d0*Mtot**2*(z+1.d0)**2*(ms**2-g)+16.d0*(g+ms**2)**2))

d11_6 = -(2.d0*Mtot*(-1.d0 + v)*(1.d0 + z)**4                                  &
        / (Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)      &
           + 16.d0*(ms**2 + g)**2))                                              &
        * ( Mtot**3*(1.d0 + z)*(-1.d0 + 3.d0*z)                                 &
            - 2.d0*Mtot**2*mf*(-1.d0 + z**2)                                    &
            - 8.d0*mf*(ms**2 + g)                                                &
            + 4.d0*Mtot*(ms**2*(1.d0 - 2.d0*z) + 3.d0*g) )

d12_6 = ((-1.d0 + v)*(1.d0 + z)**3                                              &
        / (Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)      &
           + 16.d0*(ms**2 + g)**2))                                              &
        * ( Mtot**4*(1.d0 + z)**2*((-3.d0 + z)*z + v*(1.d0 + z)**2)            &
            + 2.d0*Mtot**3*mf*(1.d0 + z)**2*(-1.d0 - z**2 + v*(1.d0 + z)**2)   &
            + 32.d0*(mf/Mtot)*(-1.d0 + v)*(ms**2 + g)**2                        &
            + 16.d0*(ms**2 + g)*(ms**2*(v - z) + (g + v*g))                     &
            + 16.d0*Mtot*mf*(1.d0 + z)*(-ms**2*(-1.d0 + v + v*z)               &
              + (v + (-1.d0 + v)*z)*g)                                           &
            + 4.d0*Mtot**2*(1.d0 + z)*(-ms**2*((-4.d0 + z)*z                    &
              + 2.d0*v*(1.d0 + z)) - 3.d0*g                                     & ! Corrigido: Sinal de +
              + 2.d0*(v + z + v*z)*g) )

d21_6 = (4.d0*Mtot**2*(-1.d0 + v)*(1.d0 + z)**4                                &
        * (Mtot*(4.d0*mf + Mtot*(-3.d0 + z))*(1.d0 + z) + 4.d0*(ms**2 + g))   &
        / (Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)      &
           + 16.d0*(ms**2 + g)**2))

d22_6 = (2.d0*(-1.d0 + v)*(1.d0 + z)**3                                        &
        / (Mtot**4*(1.d0 + z)**4 - 8.d0*Mtot**2*(1.d0 + z)**2*(ms**2 - g)      &
           + 16.d0*(ms**2 + g)**2))                                              &
        * ( -2.d0*Mtot**3*mf*(-1.d0 + z)*(1.d0 + z)**2                         & ! Corrigido: (-1.d0 + z)
            + Mtot**4*(1.d0 + z)**2*(-2.d0 + z - z**2 + v*(1.d0 + z)**2)        &
            - 8.d0*Mtot*mf*(1.d0 + z)*(ms**2 + g)                               &
            + 16.d0*(-1.d0 + v)*(ms**2 + g)**2                                  &
            + 4.d0*Mtot**2*(1.d0 + z)*(ms**2*(3.d0 - 2.d0*v*(1.d0 + z))        &
              + g + 2.d0*(v + (-1.d0 + v)*z)*g) )

    if (DEBUG_COEF) then
      if(f2   /=f2   .or.abs(f2   )>1d10)write(*,*)'[UP] f2   =',f2   ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(Du2  /=Du2  .or.abs(Du2  )>1d10)write(*,*)'[UP] Du2  =',Du2  ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_0/=d11_0.or.abs(d11_0)>1d10)write(*,*)'[UP] d11_0=',d11_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_0/=d12_0.or.abs(d12_0)>1d10)write(*,*)'[UP] d12_0=',d12_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_0/=d21_0.or.abs(d21_0)>1d10)write(*,*)'[UP] d21_0=',d21_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_0/=d22_0.or.abs(d22_0)>1d10)write(*,*)'[UP] d22_0=',d22_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_2/=d11_2.or.abs(d11_2)>1d10)write(*,*)'[UP] d11_2=',d11_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_2/=d12_2.or.abs(d12_2)>1d10)write(*,*)'[UP] d12_2=',d12_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_2/=d21_2.or.abs(d21_2)>1d10)write(*,*)'[UP] d21_2=',d21_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_2/=d22_2.or.abs(d22_2)>1d10)write(*,*)'[UP] d22_2=',d22_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_3/=d11_3.or.abs(d11_3)>1d10)write(*,*)'[UP] d11_3=',d11_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_3/=d12_3.or.abs(d12_3)>1d10)write(*,*)'[UP] d12_3=',d12_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_3/=d21_3.or.abs(d21_3)>1d10)write(*,*)'[UP] d21_3=',d21_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_3/=d22_3.or.abs(d22_3)>1d10)write(*,*)'[UP] d22_3=',d22_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_5/=d11_5.or.abs(d11_5)>1d10)write(*,*)'[UP] d11_5=',d11_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_5/=d12_5.or.abs(d12_5)>1d10)write(*,*)'[UP] d12_5=',d12_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_5/=d21_5.or.abs(d21_5)>1d10)write(*,*)'[UP] d21_5=',d21_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_5/=d22_5.or.abs(d22_5)>1d10)write(*,*)'[UP] d22_5=',d22_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_6/=d11_6.or.abs(d11_6)>1d10)write(*,*)'[UP] d11_6=',d11_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_6/=d12_6.or.abs(d12_6)>1d10)write(*,*)'[UP] d12_6=',d12_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_6/=d21_6.or.abs(d21_6)>1d10)write(*,*)'[UP] d21_6=',d21_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_6/=d22_6.or.abs(d22_6)>1d10)write(*,*)'[UP] d22_6=',d22_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
    end if

    if (DEBUG_COEF) then
    if(z>0.7d0  .and. g > 0.3d0 .and. g < 0.4d0 .and. gp >1.4d0 .and. gp<1.6d0 .and. z<0.83d0 .and. zp <0.85d0 .and. &
    zp >0.8d0 .and. v==0.5d0 .and. u == 0.25d0 .and. s ==1 .and. f==1) then
      !write(*,*) 'aaaa'
      WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] f2=', f2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] Du2  =',Du_c, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d11_0=',d11_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d12_0=',d12_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d21_0=',d21_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d22_0=',d22_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d11_2=',d11_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d12_2=',d12_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d21_2=',d21_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d22_2=',d22_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d11_3=',d11_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d12_3=',d12_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d21_3=',d21_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d22_3=',d22_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d11_5=',d11_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d12_5=',d12_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d21_5=',d21_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d22_5=',d22_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d11_6=',d11_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d12_6=',d12_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d21_6=',d21_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[UP] d22_6=',d22_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
    end if
  end if 

if (f==1 .and. s==1) Bsf_xi = (d11_0 + f2 * (d11_2 + 0.25*d11_3*(Mtot/2.0*ku - z*Mtot**2/4.0)) &
                               + 0.25*d11_5*Mtot**2 &
                               + 0.25*d11_6*(-Mtot*z/2.0*ku - g))
if (f==1 .and. s==2) Bsf_xi = (d12_0 + f2 * (d12_2 + 0.25*d12_3*(Mtot/2.0*ku - z*Mtot**2/4.0)) &
                               + 1.5*d12_4*f2 &
                               + 0.25*d12_5*Mtot**2 &
                               + 0.25*d12_6*(-Mtot*z/2.0*ku - g))

if (f==2 .and. s==1) Bsf_xi = (d21_0 + f2 * (d21_2 + 0.25*d21_3*(Mtot/2.0*ku - z*Mtot**2/4.0)) &
                               + 0.25*d21_5*Mtot**2 &
                               + 0.25*d21_6*(-Mtot*z/2.0*ku - g))

if (f==2 .and. s==2) Bsf_xi = (d22_0 + f2 * (d22_2 + 0.25*d22_3*(Mtot/2.0*ku - z*Mtot**2/4.0)) &
                               + 1.5*d22_4*f2 &
                               + 0.25*d22_5*Mtot**2 &
                               + 0.25*d22_6*(-Mtot*z/2.0*ku - g))

    denom1 = Du_c**3
    denom2 = 1.d0


    Pij2 = -2.0 * v**2 / (D0 * denom1 * denom2) * Bsf_xi

    if (DEBUG_COEF) then
      if(Bsf_xi/=Bsf_xi.or.abs(Bsf_xi)>1d10) &
        write(*,*)'[UP] Bsf_xi =',Bsf_xi,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(denom1/=denom1.or.abs(denom1)<1d-60) &
        write(*,*)'[UP] denom1~0=',denom1,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(denom2/=denom2.or.abs(denom2)<1d-60) &
        write(*,*)'[UP] denom2~0=',denom2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(Pij2  /=Pij2  .or.abs(Pij2  )>1d10) &
        write(*,*)'[UP] Pij2   =',Pij2  ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
    end if
    
    
  end if

  else
    ! Coeficientes escalares
    a     = 1.d0 + 2.d0*mf/Mtot
    b_sca = 0.5d0 - mf/Mtot
    c11_0 = 0.5d0*Mtot*a
    c12_0 = -0.25d0*zp*v*Mtot*a - (1.d0-v)/Mtot*(g + z*M2_4)
    c12_1 = 0.5d0*(1.d0-v)*(1.d0-z)
    c21_0 = Mtot
    c22_0 = -0.5d0*Mtot*zp*v - (1.d0-v)*Mtot*b_sca
    if (f==1 .and. s==1) Bsf = c11_0
    if (f==1 .and. s==2) Bsf = c12_0 + c12_1*ku
    if (f==2 .and. s==1) Bsf = c21_0
    if (f==2 .and. s==2) Bsf = c22_0
  end if
  
  
  KERNEL_UPPER = (1.d0+z)**2 / (32.d0*PI**2 * D0) * v**2 * Bsf / Du_e**2 + (1.d0-xi)* Pij2/(8.d0*PI**2)
END FUNCTION KERNEL_UPPER


! -----------------------------------------------------------------------
! KERNEL_LOWER  —  integrando para z' em [-1, z]
! -----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION KERNEL_LOWER(z, zp, g, gp, v, m, mu, kappa, PI, s, f, ms, mf, Mtot, kernel_type, xi, u)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)
  DOUBLE PRECISION, INTENT(IN) :: z, zp, g, gp, v, PI, u
  DOUBLE PRECISION, INTENT(IN) :: m, mu, kappa, ms, mf, Mtot, xi
  INTEGER, INTENT(IN) :: s, f
  CHARACTER, INTENT(IN) :: kernel_type
  DOUBLE PRECISION :: M2_4, kd, lD1, Ld2  ! k,M começam com letra i-n: precisam de declaração explícita
  LOGICAL, PARAMETER :: DEBUG_COEF = .TRUE.  ! mude para .TRUE. para ativar flags de explosão

  M2_4  = 0.25d0*Mtot**2
  Delta = 0.5d0*(ms - mf)

  ! Eq. (16): D_0(gamma, z)
  D0 = g + (1.d0 - z**2)*kappa**2 + (Delta - m*z)**2

  ! Eq. (18): D_d(z', z, mF^2) = D_u(-z', -z, mF^2)
  Dd_e = v*(1.d0-v)*(z-zp)*( g - (1.d0-z**2)*M2_4 + mf**2 )  &
     + (1.d0-z)*( v*(1.d0-v)*( g + z**2*M2_4 )              &
                + v*(gp + kappa**2)                           &
                + v**2*zp**2*M2_4                             &
                + (1.d0-v)*mu**2 )

  

  ! Eq. (20): k_d^-
  kd = -0.5d0*Mtot + 2.d0*(g + mf**2)/(Mtot*(1.d0-z))

  Bsf  = 0.d0
  Pij2 = 0.d0
  if (kernel_type == 'v') then
    ! Coeficientes vetoriais
    a     = 1.d0 + 2.d0*mf/Mtot
    c11_0 = 0.5d0*Mtot**2*(a + 0.5d0*v*zp*a + 0.5d0*(2.d0-v)*z) + (2.d0-v)*g
    c11_1 = -0.5d0*Mtot*(2.d0-v)*(1.d0-z)
    c12_0 = 0.5d0*a*( g*(1.d0-v)*(2.d0+v) + 2.d0*v*(gp+kappa**2)  &
          + 2.d0*(1.d0-v)*mu**2 + M2_4*(v*zp-2.d0)*(z-v*(z-zp)) )  &
          + (g + z*M2_4)*(1.d0-v-0.5d0*zp*v)
    c12_1 = 0.5d0*Mtot*( a*(1.d0-v)*(1.d0+z+0.5d0*(z-zp)*v)        &
          - (1.d0-z)*(1.d0-v-0.5d0*zp*v) )
    c21_0 = 0.5d0*Mtot**2*((2.d0-v)*(2.d0-a) + 2.d0 + zp*v)
    c22_0 = g*(1.d0-v)*(2.d0+v) + 2.d0*v*(gp+kappa**2) + 2.d0*(1.d0-v)*mu**2  &
          + 0.5d0*Mtot**2*( (1.d0-v)*(0.5d0*z*zp*v - z + 2.d0-a)               &
                          + 0.5d0*zp*v*(zp*v - 4.d0+a) )
    c22_1 = 0.5d0*Mtot*(1.d0-v)*(v*(z-zp) + 2.d0*(1.d0+z))
    if (f==1 .and. s==1) Bsf = c11_0 + c11_1*kd
    if (f==1 .and. s==2) Bsf = c12_0 + c12_1*kd
    if (f==2 .and. s==1) Bsf = c21_0
    if (f==2 .and. s==2) Bsf = c22_0 + c22_1*kd


  if (xi /= 1.d0) then
    !Constantes Alines

  f2 = - g * v*(1.d0 - v)*(1.d0 - zp) &
        - v*(1.d0 - z)*(gp + kappa**2) &
        + (Mtot**2 * v*(1.d0 - z) / 4.d0) &
          * (z*(1.d0 - zp)*(1.d0 - v) - zp*(1.d0 - v + v*zp)) &
        - mf**2 * v*(1.d0 - v)*(z - zp) &
        - mu**2 * (1.d0 - z)*(1.d0 - u - v + xi*u)

  Dd_calibre = (1.d0/4.d0)*v &
     * (4.d0*g + 4.d0*gp + Mtot**2*(z-1.d0)*(zp-1.d0)*((v-1.d0)*z-v*zp-1.d0) &
     + mf**2*((3.d0-4.d0*v)*z + 4.d0*(v-1.d0)*zp + 1.d0) &
     - 2.d0*mf*ms*(z-1.d0) - ms**2*z + ms**2 &
     - 4.d0*g*v + 4.d0*g*v*zp - 4.d0*gp*z - 4.d0*g*zp) &
     + mu**2*(z-1.d0)*(-xi*u+u+v-1.d0)

!Coeficientes Aline Kd

d11_0 = ((Mtot + 2.d0*mf)*v**2*(2.d0*mf + Mtot*(-1.d0 + zp))*(1.d0 - z)**2       &
        / 32.d0)                                                                  &
        * ( 4.d0*mf**2*(-1.d0 + (-2.d0 + v)*z + zp - v*zp)                     &
            - Mtot**2*(-1.d0 + z)                                                &
            * (1.d0 + z*(3.d0 - v + (-1.d0 + v)*zp) + zp*(-3.d0 + v - v*zp))   &
            + 4.d0*(-3.d0 + v + zp - v*zp)*g )

d12_0 = (v**2*(Mtot+2.d0*mf)*(1-z) &
     * (Mtot**2*(z-1.d0)*(zp-1.d0)*((v-1.d0)*z-v*zp) &
        - 2.d0*Mtot*mf*(z-1.d0)*zp &
        - 4.d0*(v-1.d0)*(g+mf**2*(z-zp)-g*zp)) &
     * (Mtot**2*(z-1.d0)*(z*(v*(zp-1.d0)-zp+3.d0)-v*zp**2+(v-3.d0)*zp+1.d0) &
        - 4.d0*(mf**2*((v-2.d0)*z-v*zp+zp-1.d0)+g*(v*(-zp)+v+zp-3.d0)))) &
     / (64.d0*Mtot)

d21_0 = (Mtot*v**2*(1.d0 - z)**2 / 16.d0)                                          &
        * ( (2.d0*mf + Mtot*(-1.d0 + zp))                                       &
            * (4.d0*mf**2*(-1.d0 + (-2.d0 + v)*z + zp - v*zp)                  &
               - Mtot**2*(-1.d0 + z)                                             &
               * (1.d0 + z*(3.d0 - v + (-1.d0 + v)*zp)                          &
                  + zp*(-3.d0 + v - v*zp))                                       &
               + 4.d0*(-3.d0 + v + zp - v*zp)*g) )

d22_0 = (v**2 / 32.d0)*(1-z)                                                  &
        * ( -4.d0*mf**2*(-1.d0 + (-2.d0 + v)*z + zp - v*zp)                    &
            + Mtot**2*(-1.d0 + z)*(1.d0 + z*(3.d0 - v + (-1.d0 + v)*zp)        &
              + zp*(-3.d0 + v - v*zp))                                           &
            + 4.d0*(3.d0 - v - zp + v*zp)*g )                                   &
        * ( -2.d0*Mtot*mf*(-1.d0 + z)*zp                                        &
            + Mtot**2*(-1.d0 + z)*(-1.d0 + zp)*((-1.d0 + v)*z - v*zp)          &
            - 4.d0*(-1.d0 + v)*(mf**2*(z - zp) + g - zp*g) )                    ! Corrigido

d11_2 = -(v*(1.d0 - z)**2 / 8.d0)                                               &
        * (Mtot + 2.d0*mf)*(2.d0*mf + Mtot*(-1.d0 + zp))

d12_2 = (v*(Mtot + 2.d0*mf)*(1.d0 - z) / (16.d0*Mtot))                         &
        * ( -2.d0*Mtot*mf*(-1.d0 + z)*zp                                        &
            + 4.d0*mf**2*(1.d0 + (3.d0 - 2.d0*v)*z + 2.d0*(-1.d0 + v)*zp)      &
            + Mtot**2*(-1.d0 + z)*(1.d0 + 2.d0*z*(2.d0 - v                      &
              + (-1.d0 + v)*zp) + zp*(-3.d0 + 2.d0*v - 2.d0*v*zp))             &
            + 8.d0*(2.d0 + v*(-1.d0 + zp) - zp)*g )

d21_2 = -(Mtot*v*(1.d0 - z)**2 / 4.d0)                                          &
        * (2.d0*mf + Mtot*(-1.d0 + zp))

d22_2 = -((1.d0 - z) / 8.d0)                                                    &
        * ( v*(2.d0*Mtot*mf*(-1.d0 + z)*zp                                      &
            + 4.d0*mf**2*(-1.d0 - 3.d0*z + 2.d0*v*z + 2.d0*zp - 2.d0*v*zp)    & ! Corrigido
            - Mtot**2*(-1.d0 + z)*(1.d0 + 2.d0*z*(2.d0 - v                      &
              + (-1.d0 + v)*zp) + zp*(-3.d0 + 2.d0*v - 2.d0*v*zp))             &
            + 8.d0*(-2.d0 + v + zp - v*zp)*g) )

  d11_3 = ((z-1.d0)**3*(Mtot+2.d0*mf)*(Mtot*(2.d0*v*(zp+3.d0)-7.d0)+2.d0*mf*(2.d0*v-3.d0))) &
     / (2.d0*Mtot**2*(z-1.d0)**2-8.d0*Mtot*mf*(z-1.d0)+8.d0*(g+mf**2))

d12_3 = ((z-1.d0)**2*(Mtot+2.d0*mf) &
     * (Mtot**2*(z-1.d0)*(4.d0*v**2*(z-1.d0)*(zp+1.d0)+v*(-6.d0*z*(zp+2.d0)+4.d0*zp+6.d0)+9.d0*z-2.d0) &
     - 4.d0*Mtot*mf*(z-1.d0)*(4.d0*v**2*(zp+1.d0)-v*(5.d0*zp+8.d0)+4.d0) &
     + 4.d0*(2.d0*v-3.d0)*(mf**2*(2.d0*v*(zp+1.d0)-z-2.d0)+g*(2.d0*v*(zp+1.d0)-3.d0)))) &
     / (4.d0*Mtot*(Mtot**2*(z-1.d0)**2-4.d0*Mtot*mf*(z-1.d0)+4.d0*(g+mf**2)))

d21_3 = (Mtot*(z-1.d0)**3*(Mtot*(2.d0*v*(zp+3.d0)-7.d0)+2.d0*mf*(2.d0*v-3.d0))) &
     / (Mtot**2*(z-1.d0)**2-4.d0*Mtot*mf*(z-1.d0)+4.d0*(g+mf**2))

d22_3 = ((z-1.d0)**2 &
     * (Mtot**2*(z-1.d0)*(4.d0*v**2*(z-1.d0)*(zp+1.d0)+v*(-6.d0*z*(zp+2.d0)+4.d0*zp+6.d0)+9.d0*z-2.d0) &
     - 4.d0*Mtot*mf*(z-1.d0)*(4.d0*v**2*(zp+1.d0)-v*(5.d0*zp+8.d0)+4.d0) &
     + 4.d0*(2.d0*v-3.d0)*(mf**2*(2.d0*v*(zp+1.d0)-z-2.d0)+g*(2.d0*v*(zp+1.d0)-3.d0)))) &
     / (2.d0*Mtot**2*(z-1.d0)**2-8.d0*Mtot*mf*(z-1.d0)+8.d0*(g+mf**2))

d12_4 = (Mtot + 2.d0*mf)*(1.d0 - z)**2                                          &
        / (4.d0*Mtot)

d22_4 = (1.d0 - z)**2 / 2.d0

d11_5 = -(((z-1.0d0)**4*(Mtot+2.0d0*mf)**2*(v*zp+1.0d0))/(2.0d0*Mtot**2*(z-1.0d0)**2 &
        -8.0d0*Mtot*mf*(z-1.0d0)+8.0d0*(g+mf**2)))

d12_5 = ((Mtot + 2.d0*mf)*(1.d0 + v*zp)*(1.d0 - z)**3                          &
        / (4.d0*((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g)))                       &
        * ( -4.d0*mf*v*(-1.d0 + z)*zp                                           &
            - 4.d0*(mf**2/Mtot)*(z - v*zp)                                      &
            + Mtot*(-1.d0 + z)*(-z + v*(-1.d0 + z)*zp)                          &
            + (4.d0/Mtot)*(-1.d0 + v*zp)*g )

d21_5 = (Mtot*(Mtot + 2.d0*mf)*(-1.d0 + z)*(1.d0 - z)**3*(1.d0 + v*zp))       &
        / ((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g)

d22_5 = ((1.d0 - z)**3                                                           &
        / (2.d0*(Mtot + 2.d0*mf - Mtot*z)**2 + 8.d0*g))                        &
        * ( (1.d0 + v*zp)*(-4.d0*Mtot*mf*v*(-1.d0 + z)*zp                      &
            - 4.d0*mf**2*(z - v*zp)                                              &
            + Mtot**2*(-1.d0 + z)*(-z + v*(-1.d0 + z)*zp)                       &
            + 4.d0*(-1.d0 + v*zp)*g) )

  d11_6 = (2.d0*Mtot*(Mtot + 2.d0*mf)*(-1.d0 + v)*(-1.d0 + z)*(1.d0 - z)**3)    &
        / ((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g)

d12_6 = ((1.d0 - z)**3                                                           &
        / (Mtot*((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g)))                       &
        * ( (Mtot + 2.d0*mf)*(-1.d0 + v) * ( Mtot**2*(v*(-1.d0 + z) - z)       & ! <-- Parêntese aberto aqui
              *(-1.d0 + z)                                                        &
            + 2.d0*Mtot*mf*(-1.d0 - 2.d0*v*(-1.d0 + z) + z)                    &
            + 4.d0*(-1.d0 + v)*(mf**2 + g) ) )                                   ! <-- Parêntese fechado aqui

d21_6 = (4.d0*Mtot**2*(-1.d0 + v)*(-1.d0 + z)*(1.d0 - z)**3)                   &
        / ((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g)

d22_6 = ((1.d0 - z)**3                                                           &
        / ((Mtot + 2.d0*mf - Mtot*z)**2 + 4.d0*g))                              &
        * ( 2.d0*(-1.d0 + v) * ( Mtot**2*(v*(-1.d0 + z) - z)*(-1.d0 + z)       & ! <-- Parêntese aberto aqui
            + 2.d0*Mtot*mf*(-1.d0 - 2.d0*v*(-1.d0 + z) + z)                    &
            + 4.d0*(-1.d0 + v)*(mf**2 + g) ) )                                   ! <-- Parêntese fechado aqui


    if (DEBUG_COEF) then
      if(f2   /=f2   .or.abs(f2   )>1d30)write(*,*)'[LO] f2   =',f2   ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      !if(Dd_calibre  /=Dd_calibre  .or.abs(Dd2  )>1d30)write(*,*)'[LO] Dd2  =',Dd_calibre  ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_0/=d11_0.or.abs(d11_0)>1d30)write(*,*)'[LO] d11_0=',d11_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_0/=d12_0.or.abs(d12_0)>1d30)write(*,*)'[LO] d12_0=',d12_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_0/=d21_0.or.abs(d21_0)>1d30)write(*,*)'[LO] d21_0=',d21_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_0/=d22_0.or.abs(d22_0)>1d30)write(*,*)'[LO] d22_0=',d22_0,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_2/=d11_2.or.abs(d11_2)>1d30)write(*,*)'[LO] d11_2=',d11_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_2/=d12_2.or.abs(d12_2)>1d30)write(*,*)'[LO] d12_2=',d12_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_2/=d21_2.or.abs(d21_2)>1d30)write(*,*)'[LO] d21_2=',d21_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_2/=d22_2.or.abs(d22_2)>1d30)write(*,*)'[LO] d22_2=',d22_2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_3/=d11_3.or.abs(d11_3)>1d30)write(*,*)'[LO] d11_3=',d11_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_3/=d12_3.or.abs(d12_3)>1d30)write(*,*)'[LO] d12_3=',d12_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_3/=d21_3.or.abs(d21_3)>1d30)write(*,*)'[LO] d21_3=',d21_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_3/=d22_3.or.abs(d22_3)>1d30)write(*,*)'[LO] d22_3=',d22_3,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_5/=d11_5.or.abs(d11_5)>1d30)write(*,*)'[LO] d11_5=',d11_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_5/=d12_5.or.abs(d12_5)>1d30)write(*,*)'[LO] d12_5=',d12_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_5/=d21_5.or.abs(d21_5)>1d30)write(*,*)'[LO] d21_5=',d21_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_5/=d22_5.or.abs(d22_5)>1d30)write(*,*)'[LO] d22_5=',d22_5,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d11_6/=d11_6.or.abs(d11_6)>1d30)write(*,*)'[LO] d11_6=',d11_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d12_6/=d12_6.or.abs(d12_6)>1d30)write(*,*)'[LO] d12_6=',d12_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d21_6/=d21_6.or.abs(d21_6)>1d30)write(*,*)'[LO] d21_6=',d21_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(d22_6/=d22_6.or.abs(d22_6)>1d30)write(*,*)'[LO] d22_6=',d22_6,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
    end if

    if (DEBUG_COEF) then
    if(z>0.7d0 .and. g > 0.3d0 .and. g < 0.4d0 .and. gp >1.4d0 .and. gp<1.6d0 .and. z<0.8d0 .and. zp <0.71d0 .and. &
    zp >0.7d0 .and. v==0.5d0 .and. u == 0.25d0 .and. s ==1 .and. f==1) then

      WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] f2=', f2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] Dd2  =',Dd_calibre, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d11_0=',d11_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d12_0=',d12_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d21_0=',d21_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d22_0=',d22_0, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d11_2=',d11_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d12_2=',d12_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d21_2=',d21_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d22_2=',d22_2, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d11_3=',d11_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d12_3=',d12_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d21_3=',d21_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d22_3=',d22_3, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d11_5=',d11_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d12_5=',d12_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d21_5=',d21_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d22_5=',d22_5, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d11_6=',d11_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d12_6=',d12_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d21_6=',d21_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
     WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
     '[LO] d22_6=',d22_6, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
    end if
    end if

if (f==1 .and. s==1) Bsf_xi = (d11_0 + f2 * (d11_2 + 0.25*d11_3*(Mtot/2.0*kd - z*Mtot**2/4.0)) &
                               + 0.25*d11_5*Mtot**2 &
                               + 0.25*d11_6*(-Mtot*z/2.0*kd - g))
if (f==1 .and. s==2) Bsf_xi = (d12_0 + f2 * (d12_2 + 0.25*d12_3*(Mtot/2.0*kd - z*Mtot**2/4.0)) &
                               + 1.5*d12_4*f2 &
                               + 0.25*d12_5*Mtot**2 &
                               + 0.25*d12_6*(-Mtot*z/2.0*kd - g))

if (f==2 .and. s==1) Bsf_xi = (d21_0 + f2 * (d21_2 + 0.25*d21_3*(Mtot/2.0*kd - z*Mtot**2/4.0)) &
                               + 0.25*d21_5*Mtot**2 &
                               + 0.25*d21_6*(-Mtot*z/2.0*kd - g))

if (f==2 .and. s==2) Bsf_xi = (d22_0 + f2 * (d22_2 + 0.25*d22_3*(Mtot/2.0*kd - z*Mtot**2/4.0)) &
                               + 1.5*d22_4*f2 &
                               + 0.25*d22_5*Mtot**2 &
                               + 0.25*d22_6*(-Mtot*z/2.0*kd - g))

    denom1 = Dd_calibre**3
    denom2 = 1.d0
    Pij2 = -2.0 * v**2 / (D0 * denom1 * denom2) * Bsf_xi

    if (DEBUG_COEF) then
      if(Bsf_xi/=Bsf_xi.or.abs(Bsf_xi)>1d30) &
        write(*,*)'[LO] Bsf_xi =',Bsf_xi,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(denom1/=denom1.or.abs(denom1)<1d-60) &
        write(*,*)'[LO] denom1~0=',denom1,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(denom2/=denom2.or.abs(denom2)<1d-60) &
        write(*,*)'[LO] denom2~0=',denom2,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
      if(Pij2  /=Pij2  .or.abs(Pij2  )>1d30) &
        write(*,*)'[LO] Pij2   =',Pij2  ,' z=',z,' zp=',zp,' g=',g,' gp=',gp,' v=',v,' s,f=',s,f
    end if
  end if

  else
    ! Coeficientes escalares
    a     = 1.d0 + 2.d0*mf/Mtot
    b_sca = 0.5d0 - mf/Mtot
    c11_0 = 0.5d0*Mtot*a
    c12_0 = -0.25d0*zp*v*Mtot*a - (1.d0-v)/Mtot*(g + z*M2_4)
    c12_1 = 0.5d0*(1.d0-v)*(1.d0-z)
    c21_0 = Mtot
    c22_0 = -0.5d0*Mtot*zp*v - (1.d0-v)*Mtot*b_sca

    if (f==1 .and. s==1) Bsf = c11_0
    if (f==1 .and. s==2) Bsf = c12_0 + c12_1*kd
    if (f==2 .and. s==1) Bsf = c21_0
    if (f==2 .and. s==2) Bsf = c22_0
  end if


  !razao = ((1.d0-xi)* Pij2) / ((1.d0-z)**2 / (32.d0*PI**2 * D0) * v**2 * Bsf / Dd**2)
!if(z>2.d0  .and. g > 0.3d0 .and. g < 0.4d0 .and. gp >1.4d0 .and. gp<1.6d0 .and. z<0.8d0 .and. zp <0.71d0 .and. &
 !   zp >0.7d0 .and. v==0.5d0 .and. u == 0.25d0 .and. s ==1 .and. f==1) then

  !    WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
   !  '[LO] razao=', razao, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
    !end if

 ! if(razao > 1d30) then
! WRITE(15,'(A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,2I6)') &
 !    '[LO] razao=', razao, ' z=', z, ' zp=', zp, ' g=', g, ' gp=', gp, ' v=', v, ' u=', u, ' s,f=', s, f
 ! end if 

KERNEL_LOWER = (1.d0-z)**2 / (32.d0*PI**2 * D0) * v**2 * Bsf / Dd_e**2 + (1.d0-xi)* Pij2/(8.d0*PI**2)

END FUNCTION KERNEL_LOWER


! -----------------------------------------------------------------------
! Rotinas de Spline
! -----------------------------------------------------------------------
      SUBROUTINE SPLGR1 (X,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),HI(NP1),U(NP1),Q(NP1,NP1),C(NP1,NP1)
      COMMON /FActz/ FAK1(NP1,NP1),FAK2(NP1,NP1),FAK3(NP1,NP1)
      U(1)=0.D0
      HI(2)=X(2)-X(1)
      DO 5 I=1,N
    5  Q(1,I)=0.D0
      DO 10 I=2,N-1
       AX=X(I+1)-X(I)
       HI(I+1)=AX
       BX=X(I+1)-X(I-1)
       CX=X(I)-X(I-1)
       AL=AX/BX
       AM=1.D0-AL
       PI=1.D0/(2.D0-AM*U(I-1))
       U(I)=AL*PI
       DO 15 J=1,N
   15   Q(I,J)=-PI*AM*Q(I-1,J)
       Q(I,I-1)=Q(I,I-1)+PI/(CX*BX)
       Q(I,I)=Q(I,I)-PI/(CX*AX)
   10  Q(I,I+1)=Q(I,I+1)+PI/(AX*BX)
      DO 20 J=1,N
       C(N,J)=0.D0
       FAK1(N,J)=0.D0
       FAK2(N,J)=0.D0
   20  FAK3(N,J)=0.D0
      DO 25 I=N-1,1,-1
       H1=1.D0/HI(I+1)
       DO 30 J=1,N
        C(I,J)=Q(I,J)-C(I+1,J)*U(I)
   30   FAK1(I,J)=-HI(I+1)*(2.D0*C(I,J)+C(I+1,J))
       FAK1(I,I)=FAK1(I,I)-H1
       FAK1(I,I+1)=FAK1(I,I+1)+H1
       DO 25 J=1,N
        FAK2(I,J)=3*C(I,J)
   25   FAK3(I,J)=(C(I+1,J)-C(I,J))*H1
      return
      END

      SUBROUTINE SPLMD1 (X,N,XA,SPL)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),SPL(N)
      COMMON /FActz/ FAK1(NP1,NP1),FAK2(NP1,NP1),FAK3(NP1,NP1)
      I=-1
   10  I=I+1
       IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10
      IF (I .EQ. 0) I=1
      DX=XA-X(I)
      DO 20 J=1,N
   20  SPL(J)=((FAK3(I,J)*DX+FAK2(I,J))*DX+FAK1(I,J))*DX
      SPL(I)=SPL(I)+1.D0
      return
      END

      SUBROUTINE SPLGR2 (X,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),HI(NP1),U(NP1),Q(NP1,NP1),C(NP1,NP1)
      COMMON /FActg/ FAK4(NP1,NP1),FAK5(NP1,NP1),FAK6(NP1,NP1)
      U(1)=0.D0
      HI(2)=X(2)-X(1)
      DO 5 I=1,N
    5  Q(1,I)=0.D0
      DO 10 I=2,N-1
       AX=X(I+1)-X(I)
       HI(I+1)=AX
       BX=X(I+1)-X(I-1)
       CX=X(I)-X(I-1)
       AL=AX/BX
       AM=1.D0-AL
       PI=1.D0/(2.D0-AM*U(I-1))
       U(I)=AL*PI
       DO 15 J=1,N
   15   Q(I,J)=-PI*AM*Q(I-1,J)
       Q(I,I-1)=Q(I,I-1)+PI/(CX*BX)
       Q(I,I)=Q(I,I)-PI/(CX*AX)
   10  Q(I,I+1)=Q(I,I+1)+PI/(AX*BX)
      DO 20 J=1,N
       C(N,J)=0.D0
       FAK4(N,J)=0.D0
       FAK5(N,J)=0.D0
   20  FAK6(N,J)=0.D0
      DO 25 I=N-1,1,-1
       H1=1.D0/HI(I+1)
       DO 30 J=1,N
        C(I,J)=Q(I,J)-C(I+1,J)*U(I)
   30   FAK4(I,J)=-HI(I+1)*(2.D0*C(I,J)+C(I+1,J))
       FAK4(I,I)=FAK4(I,I)-H1
       FAK4(I,I+1)=FAK4(I,I+1)+H1
       DO 25 J=1,N
        FAK5(I,J)=3*C(I,J)
   25   FAK6(I,J)=(C(I+1,J)-C(I,J))*H1
      return
      END

      SUBROUTINE SPLMD2 (X,N,XA,SPL)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),SPL(N)
      COMMON /FActg/ FAK4(NP1,NP1),FAK5(NP1,NP1),FAK6(NP1,NP1)
      I=-1
   10  I=I+1
       IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10
      IF (I .EQ. 0) I=1
      DX=XA-X(I)
      DO 20 J=1,N
   20  SPL(J)=((FAK6(I,J)*DX+FAK5(I,J))*DX+FAK4(I,J))*DX
      SPL(I)=SPL(I)+1.D0
      return
      END


! -----------------------------------------------------------------------
! Quadratura de Gauss-Legendre
! -----------------------------------------------------------------------
      SUBROUTINE legauss(XS,XL,N,X,DX,ZZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension X(N),DX(N)
      IF(N)10,10,20
 10   WRITE(5,600) N
 600  FORMAT(1H ,I10,' REJEITADOS PTOS.LEG-GAUSS')
      RETURN
 20   IF(N-2) 30,40,40
 30   X(1)=0.D0
      DX(1)=.5D0
      GO TO 140
 40   I=1
      G=-1.D0
      IC=(N+1)/2
 50   S=G
      T=1.D0
      U=1.D0
      V=0.D0
      DO 60 K=2,N
      A=K
      FACT1=(2.D0*A-1.D0)/A
      FACT2=(A-1.D0)/A
      P=FACT1*G*S-FACT2*T
      DP=FACT1*(S+G*U)-FACT2*V
      T=S
      S=P
      V=U
 60   U=DP
      SUM=0.D0
      IF(I-1)90,90,70
 70   IM1=I-1
      DO 80 K=1,IM1
 80   SUM=SUM+1.D0/(G-X(K))
 90   TEST=G
      G=G-P/(DP-P*SUM)
      R=DABS(TEST-G)
      IF(R.LT.ZZ)GOTO 100
      GOTO 50
 100  R=N
      X(I)=G
      DX(I)=2.D0/R/T/DP
      IF(IC-I)120,120,110
 110  FIM1=IM1
      G=G-(DP-P*SUM)/((2.D0*G*DP-A*(A+1.D0)*P)/(1.D0-G*G)-2.D0*DP*SUM-P*SUM**2+FIM1*P)
      I=I+1
      GOTO 50
 120  K0=2*IC-N+2*(N/2)+1
      IC=IC+1
      DO 130 I=IC,N
      K=K0-I
      X(I)=-X(K)
 130  DX(I)=DX(K)
 140  FACT1=(XL-XS)/2.D0
      FACT2=(XL+XS)/2.D0
      DO 150 I=1,N
      DX(I)=DX(I)*FACT1
 150  X(I)=X(I)*FACT1+FACT2
      RETURN
      END


! -----------------------------------------------------------------------
! Geração de malha 1D e pontos de colocação
! -----------------------------------------------------------------------
      SUBROUTINE G1D(IW,X0,N,A,XN,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(0:N)
      X(0)=X0
      DX=(XN-X0)/DFLOAT(N)
      IF(A.NE.1.D0)DX=(XN-X0)*(A-1.D0)/(A**N-1.D0)
      DO I=1,N
      X(I)=X(I-1)+DX
      DX=DX*A
      enddo
      X(N)=XN
      HMIN=DMIN1(X(1)-X(0),X(N)-X(N-1))
      HMAX=DMAX1(X(1)-X(0),X(N)-X(N-1))
      IF(IW.EQ.0)RETURN
      WRITE(IW,100)
  100 FORMAT(/,2X,'ONE-DOMAIN GRID (G1D) CHARACTERISED BY'&
     //,2X,4X,'Xmin',7X,'Xmax ',6X,'N ',6X,'A',6X,'Hmin',4X,'Hmax',/)
      WRITE(IW,101)X0,XN,N,A,HMIN,HMAX
  101 FORMAT(2X,D10.5,2X,D10.5,2X,I4,4X,F6.4,2X,F10.4,2X,F10.4)
      WRITE(IW,102) (X(I),I=0,N)
  102 FORMAT(/,6(2X,D15.8))
      RETURN
      END

      SUBROUTINE COLLOC(IW,NCOL,N,X,XG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:N),XG(NCOL*N)
      IF(NCOL.EQ.2)THEN
      IG=1
      DO 2 I=1,N
      A=X(I-1)
      B=X(I)
      BPA=B+A
      BMA=B-A
      U=-0.577350269189626D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.577350269189626D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
    2 CONTINUE
      ENDIF
      IF(NCOL.EQ.3)THEN
      IG=1
      DO 3 I=1,N
      A=X(I-1)
      B=X(I)
      BPA=B+A
      BMA=B-A
      U=-0.774596669241483D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.0D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.774596669241483D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
    3 CONTINUE
      ENDIF
      IF(IW.EQ.0)RETURN
      WRITE(IW,100)
  100 FORMAT(/,2X,'COLLOCATION GRID',/)
      WRITE(IW,102) (XG(I),I=1,NCOL*N)
  102 FORMAT(6(2X,D15.8))
      RETURN
      END
