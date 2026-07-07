C=============================================================
C  Portico 2D - Pdelta
C  GDL por no: (ux, uy, rz)
C  Elemento: barra de portico plano 2 nos (axial + flexao)
C
C  Entrada padrao: portico2d_pdelta.dat
C  Saida padrao:   portico2d_pdelta.out
C
C  Formato do arquivo de entrada:
C    NN, NE, NM, NBC, NLOAD
C    METHOD              (1=penalidade, 2=eliminacao)
C    IPDELTA, MAXITER, TOL  (0/1, inteiro, real)
C    Nos:        id, x, y
C    Materiais:  id, E, A, I
C    Elementos:  id, no1, no2, matid
C    CC:         no, dof(1=ux,2=uy,3=rz), valor
C    Cargas:     no, dof(1=ux,2=uy,3=rz), valor
C=============================================================
      IMPLICIT NONE
      INTEGER GDL
      INTEGER MAXN, MAXE, MAXM, MAXBC, MAXL
      PARAMETER (MAXN=200, MAXE=400, MAXM=50, MAXBC=600,
     .           MAXL=600)

      INTEGER MAX3N
      PARAMETER (MAX3N = 3*MAXN)

      INTEGER NN, NE, NM, NBC, NLOAD, METHOD
      INTEGER IPDELTA, MAXITER, ITER, CONVERGED
      INTEGER I, J, E, ID, M, NODE, MATID
      INTEGER N1, N2, NDOF, NFREE
      INTEGER G(6)

      INTEGER BCNODE(MAXBC), BCDOF(MAXBC)
      DOUBLE PRECISION BCVAL(MAXBC)
      INTEGER LDNODE(MAXL), LDDOF(MAXL)
      DOUBLE PRECISION LDVAL(MAXL)

      DOUBLE PRECISION X(MAXN), Y(MAXN)
      INTEGER E1(MAXE), E2(MAXE), EMAT(MAXE)
      DOUBLE PRECISION EA(MAXM), AA(MAXM), EI(MAXM)

      DOUBLE PRECISION KMAT(MAX3N,MAX3N), F(MAX3N), U(MAX3N)
      DOUBLE PRECISION R(MAX3N), UOLD(MAX3N)
      DOUBLE PRECISION KWORK(MAX3N,MAX3N), FWORK(MAX3N)
      DOUBLE PRECISION URED(MAX3N)
      INTEGER IS_PRESCRIBED(MAX3N), FREE_GDLS(MAX3N)

      DOUBLE PRECISION KE(6,6), KL(6,6), KGLOC(6,6), KT(6,6)
      DOUBLE PRECISION T(6,6)
      DOUBLE PRECISION UE(6), QEL(6), QTOT(6)
      DOUBLE PRECISION L, C, S, EMod, Area, Iner
      DOUBLE PRECISION NAX(MAXE), NAXNEW(MAXE)
      DOUBLE PRECISION TOL, ERR, DEN, DIFF, ABSU

      CHARACTER*120 INFILE, OUTFILE
      INTEGER OUTUNIT
C-------------------------------------------------------------

      INFILE = 'portico2d_pdelta.dat'
      OUTFILE = 'portico2d_pdelta.out'
      OUTUNIT = 20

      OPEN(10,FILE=INFILE,STATUS='OLD')
      OPEN(OUTUNIT,FILE=OUTFILE,STATUS='UNKNOWN')

C----- Leitura basica
      READ(10,*) NN, NE, NM, NBC, NLOAD
      READ(10,*) METHOD
      READ(10,*) IPDELTA, MAXITER, TOL
      NDOF = 3*NN

C----- Nos
      DO I=1,NN
         READ(10,*) NODE, X(NODE), Y(NODE)
      END DO

C----- Materiais / secoes (id, E, A, I)
      DO I=1,NM
         READ(10,*) M, EA(M), AA(M), EI(M)
      END DO

C----- Elementos (id, no1, no2, matid)
      DO E=1,NE
         READ(10,*) ID, E1(ID), E2(ID), EMAT(ID)
      END DO

C----- Condicoes de contorno
      DO I=1,NBC
         READ(10,*) BCNODE(I), BCDOF(I), BCVAL(I)
      END DO

C----- Cargas nodais
      DO I=1,NLOAD
         READ(10,*) LDNODE(I), LDDOF(I), LDVAL(I)
      END DO

      CLOSE(10)

C----- Inicializacao
      DO I=1,NDOF
         F(I)=0.0D0
         U(I)=0.0D0
         UOLD(I)=0.0D0
         R(I)=0.0D0
         IS_PRESCRIBED(I)=0
         DO J=1,NDOF
            KMAT(I,J)=0.0D0
         END DO
      END DO

      DO I=1,NBC
         IS_PRESCRIBED(GDL(BCNODE(I),BCDOF(I))) = 1
      END DO

      DO I=1,NLOAD
         F(GDL(LDNODE(I),LDDOF(I))) =
     .      F(GDL(LDNODE(I),LDDOF(I))) + LDVAL(I)
      END DO

      DO E=1,NE
         NAX(E)=0.0D0
         NAXNEW(E)=0.0D0
      END DO

C=============================================================
C  Solucao linear ou P-Delta
C=============================================================
      CONVERGED = 0
      IF (IPDELTA.EQ.0) THEN
         CALL BUILD_STIFFNESS(NN,NE,MAX3N,MAXE,MAXM,X,Y,E1,E2,EMAT,
     .        EA,AA,EI,NAX,0,KMAT)
         CALL SOLVE_WITH_BC(METHOD,NDOF,MAX3N,KMAT,F,U,NBC,
     .        BCNODE,BCDOF,BCVAL,IS_PRESCRIBED,FREE_GDLS,NFREE,
     .        KWORK,FWORK,URED)
         CALL COMPUTE_AXIALS(NE,MAXE,MAXM,X,Y,E1,E2,EMAT,EA,AA,U,
     .        NAXNEW)
         DO E=1,NE
            NAX(E)=NAXNEW(E)
         END DO
         ITER = 1
         CONVERGED = 1
      ELSE
         IF (MAXITER.LE.0) MAXITER = 30
         IF (TOL.LE.0.0D0) TOL = 1.0D-8

         DO ITER=1,MAXITER
            CALL BUILD_STIFFNESS(NN,NE,MAX3N,MAXE,MAXM,X,Y,E1,E2,
     .           EMAT,EA,AA,EI,NAX,1,KMAT)
            CALL SOLVE_WITH_BC(METHOD,NDOF,MAX3N,KMAT,F,U,NBC,
     .           BCNODE,BCDOF,BCVAL,IS_PRESCRIBED,FREE_GDLS,NFREE,
     .           KWORK,FWORK,URED)
            CALL COMPUTE_AXIALS(NE,MAXE,MAXM,X,Y,E1,E2,EMAT,EA,AA,U,
     .           NAXNEW)

            ERR = 0.0D0
            DEN = 1.0D0
            DO I=1,NDOF
               DIFF = DABS(U(I)-UOLD(I))
               ABSU = DABS(U(I))
               IF (DIFF.GT.ERR) ERR = DIFF
               IF (ABSU.GT.DEN) DEN = ABSU
            END DO
            ERR = ERR/DEN

            DO I=1,NDOF
               UOLD(I)=U(I)
            END DO
            DO E=1,NE
               NAX(E)=NAXNEW(E)
            END DO

            IF (ITER.GT.1 .AND. ERR.LE.TOL) THEN
               CONVERGED = 1
               GOTO 100
            END IF
         END DO
 100     CONTINUE
      END IF

C----- Reconstroi rigidez final e calcula reacoes
      CALL BUILD_STIFFNESS(NN,NE,MAX3N,MAXE,MAXM,X,Y,E1,E2,EMAT,
     .     EA,AA,EI,NAX,IPDELTA,KMAT)
      CALL COMPUTE_REACTIONS(NDOF,MAX3N,KMAT,F,U,R)

C=============================================================
C  Saida
C=============================================================
      WRITE(OUTUNIT,*) '========================================'
      WRITE(OUTUNIT,*) 'PORTICO 2D - RESULTADOS'
      WRITE(OUTUNIT,*) 'Arquivo de entrada: ', INFILE
      WRITE(OUTUNIT,*) 'NN, NE = ', NN, NE
      IF (METHOD.EQ.1) THEN
         WRITE(OUTUNIT,*) 'Metodo: PENALIDADE'
      ELSE
         WRITE(OUTUNIT,*) 'Metodo: ELIMINACAO DE LINHAS/COLUNAS'
      END IF
      IF (IPDELTA.EQ.1) THEN
         WRITE(OUTUNIT,*) 'Analise: P-DELTA'
         WRITE(OUTUNIT,*) 'Iteracoes: ', ITER
         WRITE(OUTUNIT,*) 'Convergiu: ', CONVERGED
      ELSE
         WRITE(OUTUNIT,*) 'Analise: LINEAR'
      END IF

      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'DESLOCAMENTOS (por no):'
      DO I=1,NN
         WRITE(OUTUNIT,'(A,I4,A,1PE12.4,A,1PE12.4,A,1PE12.4)')
     .      'No ',I,': ux=',U(GDL(I,1)),'  uy=',U(GDL(I,2)),
     .      '  rz=',U(GDL(I,3))
      END DO

      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'REACOES (apenas gdl prescritos):'
      DO I=1,NBC
         WRITE(OUTUNIT,'(A,I4,A,I2,A,1PE12.4)') 'No ',BCNODE(I),
     .      ' dof ',BCDOF(I),'  R=',R(GDL(BCNODE(I),BCDOF(I)))
      END DO

      WRITE(OUTUNIT,*) '----------------------------------------'
      WRITE(OUTUNIT,*) 'ESFORCOS LOCAIS POR ELEMENTO:'
      WRITE(OUTUNIT,*) 'Convencao: [N1,V1,M1,N2,V2,M2]'
      IF (IPDELTA.EQ.1) THEN
         WRITE(OUTUNIT,*) 'Valores exportados: TOTAL = (ke+kg)*u_local'
      ELSE
         WRITE(OUTUNIT,*) 'Valores exportados: ELASTICO = ke*u_local'
      END IF
      DO E=1,NE
         N1 = E1(E)
         N2 = E2(E)
         MATID = EMAT(E)
         EMod = EA(MATID)
         Area = AA(MATID)
         Iner = EI(MATID)
         L = DSQRT((X(N2)-X(N1))**2 + (Y(N2)-Y(N1))**2)
         C = (X(N2)-X(N1))/L
         S = (Y(N2)-Y(N1))/L

         CALL FRAME_ELASTIC_LOCAL(EMod,Area,Iner,L,KL)
         IF (IPDELTA.EQ.1) THEN
            CALL FRAME_GEOMETRIC_LOCAL(NAX(E),L,KGLOC)
         ELSE
            CALL ZERO6(KGLOC)
         END IF
         DO I=1,6
            DO J=1,6
               KT(I,J) = KL(I,J) + KGLOC(I,J)
            END DO
         END DO
         CALL TRANSFORM6(C,S,T)

         G(1)=GDL(N1,1)
         G(2)=GDL(N1,2)
         G(3)=GDL(N1,3)
         G(4)=GDL(N2,1)
         G(5)=GDL(N2,2)
         G(6)=GDL(N2,3)
         CALL LOCAL_DISP(T,U,G,UE)
         CALL MATVEC6(KL,UE,QEL)
         CALL MATVEC6(KT,UE,QTOT)

         WRITE(OUTUNIT,'(A,I4,A,1PE12.4)') 'Elem ',E,
     .      ': Nax=',NAX(E)
         IF (IPDELTA.EQ.1) THEN
            WRITE(OUTUNIT,'(6(1PE12.4,1X))') (QTOT(I),I=1,6)
         ELSE
            WRITE(OUTUNIT,'(6(1PE12.4,1X))') (QEL(I),I=1,6)
         END IF
      END DO

      CLOSE(OUTUNIT)
      STOP
      END

C=============================================================
      INTEGER FUNCTION GDL(NODE, DOF)
      IMPLICIT NONE
      INTEGER NODE, DOF
      GDL = 3*(NODE-1) + DOF
      RETURN
      END

C=============================================================
      SUBROUTINE BUILD_STIFFNESS(NN,NE,LDA,MAXE,MAXM,X,Y,E1,E2,
     .     EMAT,EA,AA,EI,NAX,USEGEO,KMAT)
      IMPLICIT NONE
      INTEGER NN,NE,LDA,MAXE,MAXM,USEGEO
      INTEGER E1(MAXE), E2(MAXE), EMAT(MAXE)
      DOUBLE PRECISION X(*), Y(*), EA(MAXM), AA(MAXM), EI(MAXM)
      DOUBLE PRECISION NAX(MAXE), KMAT(LDA,LDA)
      DOUBLE PRECISION KL(6,6), KG(6,6), KT(6,6), KE(6,6), T(6,6)
      DOUBLE PRECISION L, C, S, EMod, Area, Iner
      INTEGER I,J,E,N1,N2,MATID,G1,G2,G3,G4,G5,G6
      INTEGER GDL

      DO I=1,3*NN
         DO J=1,3*NN
            KMAT(I,J)=0.0D0
         END DO
      END DO

      DO E=1,NE
         N1=E1(E)
         N2=E2(E)
         MATID=EMAT(E)
         EMod=EA(MATID)
         Area=AA(MATID)
         Iner=EI(MATID)

         L = DSQRT((X(N2)-X(N1))**2 + (Y(N2)-Y(N1))**2)
         C = (X(N2)-X(N1))/L
         S = (Y(N2)-Y(N1))/L

         CALL FRAME_ELASTIC_LOCAL(EMod,Area,Iner,L,KL)
         IF (USEGEO.EQ.1) THEN
            CALL FRAME_GEOMETRIC_LOCAL(NAX(E),L,KG)
         ELSE
            CALL ZERO6(KG)
         END IF

         DO I=1,6
            DO J=1,6
               KT(I,J)=KL(I,J)+KG(I,J)
            END DO
         END DO

         CALL TRANSFORM6(C,S,T)
         CALL GLOBALIZE6(T,KT,KE)

         G1=GDL(N1,1)
         G2=GDL(N1,2)
         G3=GDL(N1,3)
         G4=GDL(N2,1)
         G5=GDL(N2,2)
         G6=GDL(N2,3)
         CALL ASSEMBLE6(KMAT,LDA,KE,G1,G2,G3,G4,G5,G6)
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE FRAME_ELASTIC_LOCAL(E,A,II,L,K)
      IMPLICIT NONE
      DOUBLE PRECISION E,A,II,L,K(6,6)
      DOUBLE PRECISION EA_L, EI_L2, EI_L3, EI_L
      INTEGER I,J

      DO I=1,6
         DO J=1,6
            K(I,J)=0.0D0
         END DO
      END DO

      EA_L = E*A/L
      EI_L = E*II/L
      EI_L2 = E*II/(L*L)
      EI_L3 = E*II/(L*L*L)

      K(1,1)= EA_L
      K(1,4)=-EA_L
      K(4,1)=-EA_L
      K(4,4)= EA_L

      K(2,2)= 12.0D0*EI_L3
      K(2,3)=  6.0D0*EI_L2
      K(2,5)=-12.0D0*EI_L3
      K(2,6)=  6.0D0*EI_L2

      K(3,2)=  6.0D0*EI_L2
      K(3,3)=  4.0D0*EI_L
      K(3,5)= -6.0D0*EI_L2
      K(3,6)=  2.0D0*EI_L

      K(5,2)=-12.0D0*EI_L3
      K(5,3)= -6.0D0*EI_L2
      K(5,5)= 12.0D0*EI_L3
      K(5,6)= -6.0D0*EI_L2

      K(6,2)=  6.0D0*EI_L2
      K(6,3)=  2.0D0*EI_L
      K(6,5)= -6.0D0*EI_L2
      K(6,6)=  4.0D0*EI_L
      RETURN
      END

C=============================================================
      SUBROUTINE FRAME_GEOMETRIC_LOCAL(N,L,K)
      IMPLICIT NONE
      DOUBLE PRECISION N,L,K(6,6), FAC
      INTEGER I,J

C  N positivo = tracao; N negativo = compressao.
C  Assim, compressao reduz a rigidez lateral ao somar Kg.
      DO I=1,6
         DO J=1,6
            K(I,J)=0.0D0
         END DO
      END DO

      FAC = N/(30.0D0*L)
      K(2,2)= 36.0D0*FAC
      K(2,3)=  3.0D0*L*FAC
      K(2,5)=-36.0D0*FAC
      K(2,6)=  3.0D0*L*FAC

      K(3,2)=  3.0D0*L*FAC
      K(3,3)=  4.0D0*L*L*FAC
      K(3,5)= -3.0D0*L*FAC
      K(3,6)= -1.0D0*L*L*FAC

      K(5,2)=-36.0D0*FAC
      K(5,3)= -3.0D0*L*FAC
      K(5,5)= 36.0D0*FAC
      K(5,6)= -3.0D0*L*FAC

      K(6,2)=  3.0D0*L*FAC
      K(6,3)= -1.0D0*L*L*FAC
      K(6,5)= -3.0D0*L*FAC
      K(6,6)=  4.0D0*L*L*FAC
      RETURN
      END

C=============================================================
      SUBROUTINE TRANSFORM6(C,S,T)
      IMPLICIT NONE
      DOUBLE PRECISION C,S,T(6,6)
      INTEGER I,J
      DO I=1,6
         DO J=1,6
            T(I,J)=0.0D0
         END DO
      END DO
      T(1,1)= C
      T(1,2)= S
      T(2,1)=-S
      T(2,2)= C
      T(3,3)= 1.0D0
      T(4,4)= C
      T(4,5)= S
      T(5,4)=-S
      T(5,5)= C
      T(6,6)= 1.0D0
      RETURN
      END

C=============================================================
      SUBROUTINE GLOBALIZE6(T,KL,KG)
      IMPLICIT NONE
      DOUBLE PRECISION T(6,6), KL(6,6), KG(6,6)
      INTEGER I,J,A,B
      DOUBLE PRECISION SUM

      DO I=1,6
         DO J=1,6
            SUM=0.0D0
            DO A=1,6
               DO B=1,6
                  SUM = SUM + T(A,I)*KL(A,B)*T(B,J)
               END DO
            END DO
            KG(I,J)=SUM
         END DO
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE ASSEMBLE6(KMAT,LDA,KE,G1,G2,G3,G4,G5,G6)
      IMPLICIT NONE
      INTEGER LDA,G1,G2,G3,G4,G5,G6
      DOUBLE PRECISION KMAT(LDA,*), KE(6,6)
      INTEGER I,J,G(6)

      G(1)=G1
      G(2)=G2
      G(3)=G3
      G(4)=G4
      G(5)=G5
      G(6)=G6

      DO I=1,6
         DO J=1,6
            KMAT(G(I),G(J)) = KMAT(G(I),G(J)) + KE(I,J)
         END DO
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE SOLVE_WITH_BC(METHOD,NDOF,LDA,KIN,FEXT,U,NBC,
     .     BCNODE,BCDOF,BCVAL,ISP,FREE,NFREE,KWORK,FWORK,URED)
      IMPLICIT NONE
      INTEGER METHOD,NDOF,LDA,NBC,NFREE
      INTEGER BCNODE(NBC), BCDOF(NBC), ISP(NDOF), FREE(NDOF)
      DOUBLE PRECISION BCVAL(NBC)
      DOUBLE PRECISION KIN(LDA,LDA), FEXT(NDOF), U(NDOF)
      DOUBLE PRECISION KWORK(LDA,LDA), FWORK(NDOF), URED(NDOF)
      INTEGER I,J,II,JJ,G1,G2
      INTEGER GDL
      DOUBLE PRECISION PENA

      DO I=1,NDOF
         U(I)=0.0D0
      END DO
      DO I=1,NBC
         G1=GDL(BCNODE(I),BCDOF(I))
         U(G1)=BCVAL(I)
      END DO

      IF (METHOD.EQ.1) THEN
         DO I=1,NDOF
            FWORK(I)=FEXT(I)
            DO J=1,NDOF
               KWORK(I,J)=KIN(I,J)
            END DO
         END DO

         PENA=0.0D0
         DO I=1,NDOF
            IF (DABS(KWORK(I,I)).GT.PENA) PENA=DABS(KWORK(I,I))
         END DO
         IF (PENA.EQ.0.0D0) PENA=1.0D0
         PENA=1.0D12*PENA

         DO I=1,NBC
            G1=GDL(BCNODE(I),BCDOF(I))
            KWORK(G1,G1)=KWORK(G1,G1)+PENA
            FWORK(G1)=FWORK(G1)+PENA*BCVAL(I)
         END DO
         CALL GAUSS(NDOF,LDA,KWORK,FWORK,U)

      ELSE IF (METHOD.EQ.2) THEN
         NFREE=0
         DO I=1,NDOF
            IF (ISP(I).EQ.0) THEN
               NFREE=NFREE+1
               FREE(NFREE)=I
            END IF
         END DO

         DO II=1,NFREE
            I=FREE(II)
            FWORK(II)=FEXT(I)
            DO J=1,NDOF
               IF (ISP(J).EQ.1) THEN
                  FWORK(II)=FWORK(II)-KIN(I,J)*U(J)
               END IF
            END DO
            DO JJ=1,NFREE
               J=FREE(JJ)
               KWORK(II,JJ)=KIN(I,J)
            END DO
         END DO

         CALL GAUSS(NFREE,LDA,KWORK,FWORK,URED)
         DO II=1,NFREE
            U(FREE(II))=URED(II)
         END DO
      ELSE
         WRITE(*,*) 'Erro: METHOD deve ser 1 ou 2.'
         STOP
      END IF
      RETURN
      END

C=============================================================
      SUBROUTINE COMPUTE_AXIALS(NE,MAXE,MAXM,X,Y,E1,E2,EMAT,EA,AA,
     .     U,NAX)
      IMPLICIT NONE
      INTEGER NE,MAXE,MAXM,E,N1,N2,MATID
      INTEGER E1(MAXE), E2(MAXE), EMAT(MAXE), G(6)
      DOUBLE PRECISION X(*), Y(*), EA(MAXM), AA(MAXM), U(*), NAX(MAXE)
      DOUBLE PRECISION L,C,S,T(6,6),UE(6),DU
      INTEGER GDL

      DO E=1,NE
         N1=E1(E)
         N2=E2(E)
         MATID=EMAT(E)
         L = DSQRT((X(N2)-X(N1))**2 + (Y(N2)-Y(N1))**2)
         C = (X(N2)-X(N1))/L
         S = (Y(N2)-Y(N1))/L
         CALL TRANSFORM6(C,S,T)
         G(1)=GDL(N1,1)
         G(2)=GDL(N1,2)
         G(3)=GDL(N1,3)
         G(4)=GDL(N2,1)
         G(5)=GDL(N2,2)
         G(6)=GDL(N2,3)
         CALL LOCAL_DISP(T,U,G,UE)
         DU = UE(4)-UE(1)
         NAX(E)=EA(MATID)*AA(MATID)/L*DU
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE LOCAL_DISP(T,U,G,UE)
      IMPLICIT NONE
      DOUBLE PRECISION T(6,6), U(*), UE(6), SUM
      INTEGER G(6), I,J
      DO I=1,6
         SUM=0.0D0
         DO J=1,6
            SUM=SUM+T(I,J)*U(G(J))
         END DO
         UE(I)=SUM
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE MATVEC6(A,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION A(6,6), X(6), Y(6), SUM
      INTEGER I,J
      DO I=1,6
         SUM=0.0D0
         DO J=1,6
            SUM=SUM+A(I,J)*X(J)
         END DO
         Y(I)=SUM
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE ZERO6(A)
      IMPLICIT NONE
      DOUBLE PRECISION A(6,6)
      INTEGER I,J
      DO I=1,6
         DO J=1,6
            A(I,J)=0.0D0
         END DO
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE COMPUTE_REACTIONS(NDOF,LDA,K,F,U,R)
      IMPLICIT NONE
      INTEGER NDOF,LDA,I,J
      DOUBLE PRECISION K(LDA,LDA), F(NDOF), U(NDOF), R(NDOF)
      DO I=1,NDOF
         R(I)=0.0D0
         DO J=1,NDOF
            R(I)=R(I)+K(I,J)*U(J)
         END DO
         R(I)=R(I)-F(I)
      END DO
      RETURN
      END

C=============================================================
      SUBROUTINE GAUSS(N,LDA,A,B,X)
      IMPLICIT NONE
      INTEGER N,LDA,I,J,K,IP
      DOUBLE PRECISION A(LDA,LDA), B(N), X(N)
      DOUBLE PRECISION PIV,F,SUM,AMAX,TMP

C----- Eliminacao de Gauss com pivotamento parcial simples
      DO K=1,N-1
         IP=K
         AMAX=DABS(A(K,K))
         DO I=K+1,N
            IF (DABS(A(I,K)).GT.AMAX) THEN
               AMAX=DABS(A(I,K))
               IP=I
            END IF
         END DO

         IF (AMAX.LT.1.0D-30) THEN
            WRITE(*,*) 'Pivo ~ 0 em K=',K,'; sistema singular.'
            STOP
         END IF

         IF (IP.NE.K) THEN
            DO J=K,N
               TMP=A(K,J)
               A(K,J)=A(IP,J)
               A(IP,J)=TMP
            END DO
            TMP=B(K)
            B(K)=B(IP)
            B(IP)=TMP
         END IF

         PIV=A(K,K)
         DO I=K+1,N
            F=A(I,K)/PIV
            A(I,K)=0.0D0
            DO J=K+1,N
               A(I,J)=A(I,J)-F*A(K,J)
            END DO
            B(I)=B(I)-F*B(K)
         END DO
      END DO

      IF (DABS(A(N,N)).LT.1.0D-30) THEN
         WRITE(*,*) 'Pivo final ~ 0; sistema singular.'
         STOP
      END IF

      X(N)=B(N)/A(N,N)
      DO I=N-1,1,-1
         SUM=0.0D0
         DO J=I+1,N
            SUM=SUM+A(I,J)*X(J)
         END DO
         X(I)=(B(I)-SUM)/A(I,I)
      END DO
      RETURN
      END
