IMPLICIT DOUBLE PRECISION(a-h,o-z)

    INTEGER :: NMG, NMZ, NMA, LWORK
    DOUBLE PRECISION, ALLOCATABLE :: XMATRIX(:,:), ZMATRIX(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: VR(:,:), VL(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:), WORK(:)
    DOUBLE PRECISION, ALLOCATABLE :: zv(:), splz(:), splg(:), gv(:)
    DOUBLE PRECISION, ALLOCATABLE :: dzv(:), dgv(:), ALPHAR(:), ALPHAI(:)
    DOUBLE PRECISION, ALLOCATABLE :: BETA(:), c(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: XG(:), YG(:)

    INTEGER, ALLOCATABLE :: ipvt(:)
    COMMON/PARAM/PI
    COMMON/ALPHAINT/X(1000),DX(1000), Y(1000),DY(1000) &
    , W(1000),DW(1000), Nz, Ng, Nv, det (1000)
    INTEGER :: ii, jj, i, j, k, l, s &
    ,index1, index2, p, q, r, NPARAM, N_intervalZ, NCOL, N_intervalG
    INTEGER :: INFO, LDVL, IPRINTEIGEN, N_PLOT, k_max, k_min
    DOUBLE PRECISION :: LAMBDAR,  LAMBDAI, e, m, Mtot, mu, kappa, gam0, soma &
        , IntegralV, num2, max_gamma_visualizacao

    INTEGER :: Nnz (100), Nng (100), Nnv (100)

    DOUBLE PRECISION, EXTERNAL :: f_map, Jacobian_map, inverse_map


    open (unit = 10, file = "autovalores.dat",STATUS="UNKNOWN")
    open (unit = 11, file = "g_gamma.dat",STATUS="UNKNOWN")     ! g1,g2 vs gamma
    open (unit = 15, file = "g_z.dat",STATUS="UNKNOWN")         ! g1,g2 vs z
    open (unit = 12, file = "alfa.dat",STATUS="UNKNOWN")
    open (unit = 13, file = "plotalfa.dat",STATUS="UNKNOWN")
    open (unit = 14, file = 'erros.dat', status='unknown')
    open (UNIT = 20, FILE = "inputs.dat", STATUS="UNKNOWN")

    e = 0.0001d0
    PI = DACOS(-1.D0)       !3.14159265358979323846264338

        !Parâmetros
        !Massas
        Mtot = 1.0d0
        m = 1.0d0
        mu = 0.50d0
        kappa = sqrt(m**2 - 0.25*Mtot**2)
        xi = 1.0d0

        gam0 = 3.0d0

        Nz = 60
        Ng = 60
        Nv = 60

        READ(20,*) NPARAM
        DO I = 1, NPARAM
            READ(20,*) Nng(i), Nnz(i)
        END DO
        CLOSE(20)


    do ii = 1, NPARAM

        NMG = Nng(ii)
        NMZ = Nnz(ii)

        NMA = NMG * NMZ
        LWORK = 10 * NMA

        !ALOCAR VARIÁVEIS
        ALLOCATE( XMATRIX(NMA, NMA), ZMATRIX(NMA, NMA))
        ALLOCATE( VR(NMA, NMA), VL(NMA, 2*NMA) )
        ALLOCATE( WR(NMA), WI(NMA), WORK(LWORK) )
        ALLOCATE( zv(NMZ + 1), gv(NMG+1) )
        ALLOCATE( splz(NMZ), splg(NMG))
        ALLOCATE( ALPHAR(NMA), ALPHAI(NMA), BETA(NMA) )
        ALLOCATE( c(2, NMA, NMA) )
        ALLOCATE( XG(NMA), YG(NMA) )

        ALLOCATE( ipvt(NMA) )

    ! Ler matriz de coeficientes (unidade 16 só para leitura)
        open (unit = 16, file = "coeficientes.dat",STATUS="OLD")
        DO s = 1, 2
            DO i = 1, NMG
                READ(16,*) (c(s,i,j), j=1, NMZ)
            END DO
        END DO
        CLOSE(16)

        iw = 14
        N_intervalZ = (NMZ-1)/2
        N_intervalG = (NMG-1)/2
        NCOL = 2

        !Contrução das malhas

            call G1D(IW,-1.d0, N_intervalZ, 1.d0, 1.d0, X)
            call COLLOC(IW,2,N_intervalZ,X,XG)

            do i = 1, 2*N_intervalZ
            zv(i+1) = XG(i)
            end do

            zv(1)   = -0.999999999d0
            zv(NMZ) =  0.99999999d0

            call G1D(IW,0.d0, N_intervalG, 1.0d0, 3.d0, Y)
            call COLLOC(IW,2,N_intervalG,Y,YG)

            do i=1, 2*N_intervalG
            gv(i+1)=YG(i)
            enddo

            gv(1)   = 0.d0
            gv(NMG) = gam0

    if (xi /= 1.d0) gv(1) = 0.001d0


            call SPLGR1 (zv,Nmz)
            call SPLGR2 (gv,Nmg)


    g1_00 = 0.d0
    call SPLMD1(zv, Nmz, 0.d0, splz)
    call SPLMD2(gv, Nmg, 0.d0, splg)

      do j = 1, Nmz
          do i = 1, Nmg
              g1_00 = g1_00 + c(1,i,j) * splg(i) * splz(j)
          end do
      end do

      !=========================================================
      ! PLOT de g1 e g2 em função de GAMMA, para z fixo
      !   unit 11 -> colunas: gamma  g1  g2
      !=========================================================
      z_fixo = 0.d0
      max_gamma_visualizacao = 3.0d0
      N_PLOT = 1000
      call SPLMD1(zv, Nmz, z_fixo, splz)

          do p = 0, N_PLOT
              gamma_plot = (dble(p) / dble(N_PLOT)) * max_gamma_visualizacao

              call SPLMD2(gv, Nmg, gamma_plot, splg)

              soma1 = 0.d0
              soma2 = 0.d0
              do j = 1, Nmz
                  do i = 1, Nmg
                      soma1 = soma1 + c(1,i,j) * splg(i) * splz(j)
                      soma2 = soma2 + c(2,i,j) * splg(i) * splz(j)
                  end do
              end do

             write(11, '(ES25.17E3,2X,ES25.17E3,2X,ES25.17E3)') gamma_plot, soma1 / g1_00, soma2 / g1_00
          end do

      !=========================================================
      ! PLOT de g1 e g2 em função de Z, para gamma fixo
      !   unit 15 -> colunas: z  g1  g2
      !=========================================================
      gamma_fixo = 0.d0
      N_PLOT = 1500
      call SPLMD2(gv, Nmg, gamma_fixo, splg)

      do p = 0, N_PLOT
          z_plot = -cos( (dble(p) / dble(N_PLOT)) * PI )

          call SPLMD1(zv, Nmz, z_plot, splz)

          soma1 = 0.d0
          soma2 = 0.d0
          do j = 1, Nmz
              do i = 1, Nmg
                  soma1 = soma1 + c(1,i,j) * splg(i) * splz(j)
                  soma2 = soma2 + c(2,i,j) * splg(i) * splz(j)
              end do
          end do

          write(15, '(ES25.17E3,2X,ES25.17E3,2X,ES25.17E3)') z_plot, soma1/g1_00, soma2/g1_00
      end do

      ! =============================================================
      ! RECONSTRUÇÃO DA AMPLITUDE psi_i(gamma, xi)   -- Eq. (22)
      !
      !   psi_i(gamma,z) ~  int_0^inf dgamma'  g_i(gamma',z)
      !                     / [gamma + gamma' + m^2 z^2 + (1-z^2)kappa^2]^2
      !
      !   Normaliza-se dividindo por psi_i no ponto de referência fixo
      !   (gamma_ref, z_ref) para cancelar o fator global -i/M.
      !   Relação: xi = (1-z)/2  =>  z = 1 - 2*xi
      ! =============================================================

      ! Ponto de referência para a normalização (mesmo para g1 e g2)
      gamma_ref = 0.0d0
      xi_ref    = 0.5d0
      z_ref     = 1.d0 - 2.d0*xi_ref     ! = 0.0

      ! Pontos de Gauss em gamma' no intervalo [0, 3]
      CALL legauss(0.d0, 3.d0, Ng, Y, dY, 1.d-15)

      ! ---- Denominador da normalização: psi1_ref e psi2_ref ----
      call SPLMD1(zv, Nmz, z_ref, splz)

      psi1_ref = 0.d0
      psi2_ref = 0.d0
      do p = 1, Ng
          gammap = Y(p)
          dgp    = dY(p)

          call SPLMD2(gv, Nmg, gammap, splg)

          g1_val = 0.d0
          g2_val = 0.d0
          do j = 1, Nmz
              do i = 1, Nmg
                  g1_val = g1_val + c(1,i,j) * splg(i) * splz(j)
                  g2_val = g2_val + c(2,i,j) * splg(i) * splz(j)
              end do
          end do

          D_ref = gamma_ref + gammap + (m**2)*(z_ref**2) &
                  + (1.d0 - z_ref**2)*(kappa**2)

          psi1_ref = psi1_ref + (g1_val * dgp) / (D_ref**2)
          psi2_ref = psi2_ref + (g2_val * dgp) / (D_ref**2)
      end do

      ! =============================================================
      ! PLOT de psi1 e psi2 NORMALIZADAS em função de GAMMA (xi fixo)
      !   unit 18 -> colunas: gamma  psi1_norm  psi2_norm
      ! =============================================================
      open(unit = 18, file = "plot_psi_gamma.dat", STATUS="UNKNOWN")

      xi_num = 0.2d0
      z_num  = 1.d0 - 2.d0*xi_num
      max_gamma_visualizacao = 3.d0
      N_PLOT = 1000

      call SPLMD1(zv, Nmz, z_num, splz)

      do k_plot = 0, N_PLOT
          gamma_plot = (dble(k_plot) / dble(N_PLOT)) * max_gamma_visualizacao

          psi1_num = 0.d0
          psi2_num = 0.d0
          do p = 1, Ng
              gammap = Y(p)
              dgp    = dY(p)

              call SPLMD2(gv, Nmg, gammap, splg)

              g1_val = 0.d0
              g2_val = 0.d0
              do j = 1, Nmz
                  do i = 1, Nmg
                      g1_val = g1_val + c(1,i,j) * splg(i) * splz(j)
                      g2_val = g2_val + c(2,i,j) * splg(i) * splz(j)
                  end do
              end do

              D_num = gamma_plot + gammap + (m**2)*(z_num**2) &
                      + (1.d0 - z_num**2)*(kappa**2)

              psi1_num = psi1_num + (g1_val * dgp) / (D_num**2)
              psi2_num = psi2_num + (g2_val * dgp) / (D_num**2)
          end do

          psi1_norm = psi1_num / psi1_ref
          psi2_norm = psi2_num / psi2_ref

          write(18, '(ES25.17E3,2X,ES25.17E3,2X,ES25.17E3)') gamma_plot, psi1_norm, psi2_norm
      end do

      close(18)

      ! =============================================================
      ! PLOT de psi1 e psi2 NORMALIZADAS em função de Z (gamma fixo)
      !   unit 19 -> colunas: z  psi1_norm  psi2_norm
      ! =============================================================
      open(unit = 19, file = "plot_psi_z.dat", STATUS="UNKNOWN")

      gamma_num = 0.d0
      N_PLOT = 1500

      do k_plot = 0, N_PLOT
          z_plot = -cos( (dble(k_plot) / dble(N_PLOT)) * PI )

          call SPLMD1(zv, Nmz, z_plot, splz)

          psi1_num = 0.d0
          psi2_num = 0.d0
          do p = 1, Ng
              gammap = Y(p)
              dgp    = dY(p)

              call SPLMD2(gv, Nmg, gammap, splg)

              g1_val = 0.d0
              g2_val = 0.d0
              do j = 1, Nmz
                  do i = 1, Nmg
                      g1_val = g1_val + c(1,i,j) * splg(i) * splz(j)
                      g2_val = g2_val + c(2,i,j) * splg(i) * splz(j)
                  end do
              end do

              D_num = gamma_num + gammap + (m**2)*(z_plot**2) &
                      + (1.d0 - z_plot**2)*(kappa**2)

              psi1_num = psi1_num + (g1_val * dgp) / (D_num**2)
              psi2_num = psi2_num + (g2_val * dgp) / (D_num**2)
          end do

          psi1_norm = psi1_num / psi1_ref
          psi2_norm = psi2_num / psi2_ref

          write(19, '(ES25.17E3,2X,ES25.17E3,2X,ES25.17E3)') z_plot, psi1_norm, psi2_norm
      end do

      close(19)


      DEALLOCATE(XMATRIX, ZMATRIX, WORK, VR, VL, WR, WI, splz)
      DEALLOCATE(splg, ALPHAR, ALPHAI, BETA, c, XG, YG, gv, zv)
      DEALLOCATE(IPVT)


      end do

      close (15)
      close (11)
      CLOSE(10)
      CLOSE(12)
      close(13)
10     FORMAT(11E12.4)
18     format(5e15.6)
20     FORMAT(A70)


      close (14)
       Close(2)
    END



      !Rotinas
        !Spline 1
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

        !Spline 2
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









      SUBROUTINE legauss(XS,XL,N,X,DX,ZZ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension X(N),DX(N)
      IF(N)10,10,20
 10   WRITE(5,600) N
      WRITE(2,600) N
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


        !=====================================================================
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
!=====================================================================
      SUBROUTINE COLLOC(IW,NCOL,N,X,XG)
!
!     RETURN ABCISSES OF NCOL=2,3 GAUSS COLLOCATION POINTS
!     ON EACH OF THE N INTERVALS OF THE GRID X(0),X(1),...,X(N)
!
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