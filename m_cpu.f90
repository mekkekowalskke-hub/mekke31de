MODULE M_CPU

	DOUBLE PRECISION, PARAMETER :: PI=3.1415926535897932384626433832795D0
	
	CONTAINS
	
	SUBROUTINE INTEGRAL(SD,OD,D1,D2,D3,NUMDR,NUMKP,NUMEL,HFKS,NEN,IND,LDGG, &
	                    NUMEIG, &
	                    IW,JAC,OCC,PSI,KWEI,X,DPHI_G, &
	                    RHO,NUM0,CHECK,THECK,PHECK)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD,OD
		INTEGER, INTENT(IN) :: D1,D2,D3
		INTEGER, INTENT(IN) :: NUMDR
		INTEGER, INTENT(IN) :: NUMKP
		INTEGER, INTENT(IN) :: NUMEL
		INTEGER, INTENT(IN) :: HFKS
		INTEGER, INTENT(IN) :: NEN
		INTEGER, INTENT(IN) :: IND
		INTEGER, INTENT(IN) :: LDGG
		
		INTEGER, INTENT(IN) :: NUMEIG(HFKS)
		
		DOUBLE PRECISION, INTENT(IN) :: IW(D1)
		DOUBLE PRECISION, INTENT(IN) :: JAC(SD)
		DOUBLE PRECISION, INTENT(IN) :: OCC(SUM(NUMEIG),NUMKP)
		DOUBLE COMPLEX,   INTENT(IN) :: PSI(NEN,SUM(NUMEIG),NUMKP)
		DOUBLE PRECISION, INTENT(IN) :: KWEI(NUMKP)
		DOUBLE PRECISION, INTENT(IN) :: X(NEN,SD)
		DOUBLE PRECISION, INTENT(IN) :: DPHI_G(NUMEL*D1*D2*D3*NEN,0:SD)
		
		DOUBLE PRECISION, INTENT(OUT) :: RHO(NUMDR,D1*D2*D3)
		DOUBLE PRECISION, INTENT(OUT) :: NUM0(HFKS)
		
		DOUBLE PRECISION, INTENT(OUT) :: CHECK(NUMEL*D1*D2*D3*NEN,0:SD)
		DOUBLE PRECISION, INTENT(OUT) :: THECK(NUMEL*D1*D2*D3*SD,SD)
		DOUBLE COMPLEX,   INTENT(OUT) :: PHECK(NUMEL*D1*D2*D3*SUM(NUMEIG),NUMKP,0:SD)
		
		DOUBLE PRECISION :: DPHI(NEN,SD)
		DOUBLE PRECISION :: PHI(NEN)
		DOUBLE PRECISION :: GRADPHI(NEN,SD)
		DOUBLE PRECISION :: INVFXT(SD,SD)
		DOUBLE PRECISION :: INVFX(SD,SD)
		DOUBLE COMPLEX   :: PSI0(SUM(NUMEIG),NUMKP)
		DOUBLE COMPLEX   :: DPSI0(SUM(NUMEIG),NUMKP,SD)
		DOUBLE PRECISION :: Q1(SD),Q2(SD)
		DOUBLE PRECISION :: RHO0(NUMDR)
		
		DOUBLE PRECISION :: DETFX,DETJ
		DOUBLE PRECISION :: FX(SD,SD)
		
		INTEGER :: I,J,K,L,P
		INTEGER :: K0,K1,K2,K3
		INTEGER :: IB
		INTEGER :: POS
		
		POS = (IND-1)*D1*D2*D3*NEN
		
		IB = 0
		DO K1 = 1, D1
		  DO K2 = 1, D2
		    DO K3 = 1, D3
		      
		      IB = IB + 1
		      
		      PHI(:) = DPHI_G(POS+(IB-1)*NEN+1:POS+IB*NEN,0)
		      DPHI(:,:) = DPHI_G(POS+(IB-1)*NEN+1:POS+IB*NEN,1:SD)
		      
		      DO J = 1,NEN
		        DPHI(J,:) = DPHI(J,:) * JAC(:)
		      END DO
		      
		      CALL GRAD(SD,NEN,X,DPHI,FX)
		      
		      IF (LDGG .EQ. 2) THEN
		        CALL INV(SD,FX,INVFX,DETFX)
		        CALL MAT_TRANS(INVFX,SD,SD,INVFXT)
		        DO J = 1,NEN
		          Q1(:) = DPHI(J,:)
		          CALL MATVEC_PRO(INVFXT,SD,SD,Q1,Q2)
		          GRADPHI(J,:) = Q2(:)
		        END DO
		      END IF
		      
		      PSI0(:,:) = DCMPLX(0.0D0)
		      DO I = 1,NEN
		        PSI0(:,:) = PSI0(:,:) + PHI(I) * PSI(I,:,:)
		      END DO
		      
		      IF (LDGG .EQ. 2) THEN
		        DPSI0(:,:,:) = DCMPLX(0.0D0)
		        DO I = 1,NEN
		          DO K = 1,SD
		            DPSI0(:,:,K) = DPSI0(:,:,K) + GRADPHI(I,K) * PSI(I,:,:)
		          END DO
		        END DO
		      END IF
		      
		      RHO0(:) = 0.0D0
		      DO L = 1, HFKS
		        K0 = 0
		        IF (L .EQ. 2) K0 = NUMEIG(1)
		        DO K = K0+1, K0+NUMEIG(L)
		          DO P = 1, NUMKP
		            RHO0(L) = RHO0(L) + KWEI(P) * OCC(K,P) * ABS(PSI0(K,P))**2.0D0
		          END DO
		        END DO
		        IF (LDGG .EQ. 2) THEN
		          DO K = K0+1,K0+NUMEIG(L)
		            DO P = 1, NUMKP
		              RHO0(HFKS+(L-1)*SD+1:HFKS+L*SD) = &
		              RHO0(HFKS+(L-1)*SD+1:HFKS+L*SD) + &
		              KWEI(P) * OCC(K,P) * 2.0D0 * &
		              (REAL (PSI0(K,P)) * REAL (DPSI0(K,P,:)) + &
		               AIMAG(PSI0(K,P)) * AIMAG(DPSI0(K,P,:)))
		            END DO
		          END DO
		        END IF
		      END DO
		      
		      RHO(:,(K1-1)*D2*D3+(K2-1)*D3+K3) = RHO0(:)
		      
		      DETJ = IW(K1)*DETFX
		      IF (OD .EQ. 2) DETJ = IW(K2)*DETJ
		      DETJ = IW(K3)*DETJ
		      NUM0(:) = NUM0(:) + DETJ * RHO0(1:HFKS)
		      
		      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		      
		      PHECK((IND-1)*D1*D2*D3*SUM(NUMEIG)+(IB-1)*SUM(NUMEIG)+1: &
		         (IND-1)*D1*D2*D3*SUM(NUMEIG)+IB*SUM(NUMEIG),1,0) = PSI0(:,1)
		      DO I = 1,SD
		        PHECK((IND-1)*D1*D2*D3*SUM(NUMEIG)+(IB-1)*SUM(NUMEIG)+1: &
		           (IND-1)*D1*D2*D3*SUM(NUMEIG)+IB*SUM(NUMEIG),1,I) = DPSI0(:,1,I)
		      END DO
		      
		      THECK((IND-1)*D1*D2*D3*SD+(IB-1)*SD+1:(IND-1)*D1*D2*D3*SD+IB*SD,:)=INVFX(:,:)
		      
		      CHECK(POS+(IB-1)*NEN+1:POS+IB*NEN,0) = PHI(:)
		      CHECK(POS+(IB-1)*NEN+1:POS+IB*NEN,1:SD) = DPHI(:,:)
		      
		    END DO
		  END DO
		END DO
		
	END SUBROUTINE INTEGRAL
	
	SUBROUTINE BASIS(SD,OD,D1,D2,D3,DIM_CP,DIM_KP,DIM_OP,NEN,NUMEL,IND,D123, &
	                    SPAN,CWEI, &
	                    DPHI_G,DERIV, &
	                    TEMP,WP,BSS_BASIS)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD,OD
		INTEGER, INTENT(IN) :: D1,D2,D3
		INTEGER, INTENT(IN) :: DIM_CP(SD),DIM_KP(SD),DIM_OP(SD)
		INTEGER, INTENT(IN) :: NEN
		INTEGER, INTENT(IN) :: NUMEL
		INTEGER, INTENT(IN) :: IND
		INTEGER, INTENT(IN) :: D123
		
		INTEGER, INTENT(IN) :: SPAN(NUMEL,SD)
		
		DOUBLE PRECISION, INTENT(IN) :: CWEI(0:DIM_CP(1),0:DIM_CP(OD),0:DIM_CP(SD))
		
		DOUBLE PRECISION, INTENT(INOUT) :: DPHI_G(NUMEL*D1*D2*D3*NEN,0:SD)
		DOUBLE PRECISION, INTENT(INOUT) :: DERIV(NUMEL*SD)
		
		DOUBLE PRECISION, INTENT(OUT) :: TEMP(NUMEL,0:SD)
		DOUBLE PRECISION, INTENT(OUT) :: WP(NUMEL,0:SD)
		DOUBLE PRECISION, INTENT(IN) :: BSS_BASIS(NUMEL*D123, &
		   2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+2*(DIM_OP(3)+1))
		
		INTEGER :: K1,K2,K3
		INTEGER :: IB
		
		IB = 0
		
		DO K1 = 1, D1
		  DO K2 = 1, D2
		    DO K3 = 1, D3
		      
		      IB = IB+1
		      
		      CALL BASIS_IN(SD,OD,DIM_CP,DIM_KP,DIM_OP,NEN,NUMEL,IND,IB,D123, &
	                         SPAN,CWEI, &
	                         DPHI_G,DERIV, &
	                         TEMP,WP,BSS_BASIS)
		      
		    END DO
		  END DO
		END DO
		
	END SUBROUTINE BASIS
	
	SUBROUTINE BASIS_IN(SD,OD,DIM_CP,DIM_KP,DIM_OP,NEN,NUMEL,IND,IB,D123, &
	                 SPAN,CWEI, &
	                 DPHI_G,DERIV, &
	                 TEMP,WP,BSS_BASIS)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD,OD
		INTEGER, INTENT(IN) :: NEN
		INTEGER, INTENT(IN) :: NUMEL
		INTEGER, INTENT(IN) :: IND	! INDICATES IN WHICH INDEX WE ARE AMONG THE NUMELs
		INTEGER, INTENT(IN) :: IB,D123
		
		INTEGER, INTENT(IN) :: DIM_CP(SD),DIM_KP(SD),DIM_OP(SD)
		INTEGER, INTENT(IN) :: SPAN(NUMEL,SD)
		
		DOUBLE PRECISION, INTENT(IN) :: CWEI(0:DIM_CP(1),0:DIM_CP(OD),0:DIM_CP(SD))
		
		DOUBLE PRECISION, INTENT(INOUT) :: DPHI_G(NUMEL*D123*NEN,0:SD)
		DOUBLE PRECISION, INTENT(INOUT) :: DERIV(NUMEL*SD)
		
		DOUBLE PRECISION, INTENT(OUT) :: TEMP(NUMEL,0:SD)
		DOUBLE PRECISION, INTENT(OUT) :: WP(NUMEL,0:SD)
		DOUBLE PRECISION, INTENT(IN) :: BSS_BASIS(NUMEL*D123, &
		   2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+2*(DIM_OP(3)+1))
		
		INTEGER :: NO
		INTEGER :: J1,J2,J3
		INTEGER :: K1,K2,K3
		INTEGER :: I,J
		
		INTEGER :: BSS_IND
		
		BSS_IND = (IND-1)*D123 + IB
		
		DO NO = 1,SD
		  DERIV(SD*(IND-1)+NO) = DBLE(MOD(IND,SD))*DBLE(NO)
		END DO
		
		DO J = 0,SD
		  WP(IND,J) = 0.0D0
		END DO
		
		DO J1 = 0,DIM_OP(1)
		  K1 = SPAN(IND,1) - DIM_OP(1) + J1
		  DO J2 = 0,DIM_OP(2)
		    K2 = 0
		    IF (OD .EQ. 2) K2 = SPAN(IND,OD) - DIM_OP(2) + J2
		    DO J3 = 0,DIM_OP(3)
		      K3 = 0
		      K3 = SPAN(IND,SD) - DIM_OP(3) + J3
		      
		      TEMP(IND,0) = BSS_BASIS(BSS_IND,(J1+1))
		      IF (OD .EQ. 2) TEMP(IND,0) = TEMP(IND,0) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1))
		      TEMP(IND,0) = TEMP(IND,0) * BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(J3+1))
		      
		      TEMP(IND,1) = BSS_BASIS(BSS_IND,(DIM_OP(1)+1)+(J1+1))
		      IF (OD .EQ. 2) THEN
		        TEMP(IND,1) = TEMP(IND,1) * BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1))
		        TEMP(IND,OD) = BSS_BASIS(BSS_IND,(J1+1)) * &
		           BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(DIM_OP(2)+1)+(J2+1))
		      END IF
		      
		      TEMP(IND,1:OD) = TEMP(IND,1:OD) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(J3+1))
		      TEMP(IND,SD) = &
		         BSS_BASIS(BSS_IND,(J1+1)) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1)) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(DIM_OP(3)+1)+(J3+1))
		      
		      TEMP(IND,:) = TEMP(IND,:) * CWEI(K1,K2,K3)
		      
		      WP(IND,:) = WP(IND,:) + TEMP(IND,:)
		      
		    END DO
		  END DO
		END DO
		
		NO = 0
		DO J1 = 0,DIM_OP(1)
		  K1 = SPAN(IND,1) - DIM_OP(1) + J1
		  DO J2 = 0,DIM_OP(2)
		    K2 = 0
		    IF (OD .EQ. 2) K2 = SPAN(IND,OD) - DIM_OP(2) + J2
		    DO J3 = 0,DIM_OP(3)
		      K3 = 0
		      K3 = SPAN(IND,SD) - DIM_OP(3) + J3
		      
		      NO = NO + 1
		
		      TEMP(IND,0) = BSS_BASIS(BSS_IND,(J1+1))
		      IF (OD .EQ. 2) TEMP(IND,0) = TEMP(IND,0) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1))
		      TEMP(IND,0) = TEMP(IND,0) * BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(J3+1))
		      
		      TEMP(IND,1) = BSS_BASIS(BSS_IND,(DIM_OP(1)+1)+(J1+1))
		      IF (OD .EQ. 2) THEN
		        TEMP(IND,1) = TEMP(IND,1) * BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1))
		        TEMP(IND,OD) = BSS_BASIS(BSS_IND,(J1+1)) * &
		           BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(DIM_OP(2)+1)+(J2+1))
		      END IF
		      
		      TEMP(IND,1:OD) = TEMP(IND,1:OD) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(J3+1))
		      TEMP(IND,SD) = &
		         BSS_BASIS(BSS_IND,(J1+1)) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+(J2+1)) * &
		         BSS_BASIS(BSS_IND,2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+(DIM_OP(3)+1)+(J3+1))
		      
		      TEMP(IND,:) = TEMP(IND,:) * CWEI(K1,K2,K3)
		      
		      DPHI_G((IND-1)*D123*NEN+(IB-1)*NEN+NO,0) = TEMP(IND,0) / WP(IND,0)
		      DPHI_G((IND-1)*D123*NEN+(IB-1)*NEN+NO,1:SD) = &
		      ( TEMP(IND,1:SD) - WP(IND,1:SD)*DPHI_G((IND-1)*D123*NEN+(IB-1)*NEN+NO,0) ) / &
		      WP(IND,0)
		      
		    END DO
		  END DO
		END DO
		
		DO NO = 1,SD
		  DPHI_G((IND-1)*D123*NEN+(IB-1)*NEN+1:(IND-1)*D123*NEN+(IB-1)*NEN+NEN,NO) = &
		     DPHI_G((IND-1)*D123*NEN+(IB-1)*NEN+1:(IND-1)*D123*NEN+(IB-1)*NEN+NEN,NO) * &
		     DERIV(SD*(IND-1)+NO)
		END DO
		
	END SUBROUTINE BASIS_IN
	
	SUBROUTINE DERS_BASIS_FUNS(SD,OD,NUMEL,IND,IB,D123,N,MOP, &
	                           DIM_KP,SPAN,BAK_KNOT0,DIM_OP,BAK_KNOT,BSS_BASIS, &
	                           BSS_A,BSS_NDU,LEFT,RIGHT)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD,OD
		INTEGER, INTENT(IN) :: NUMEL
		INTEGER, INTENT(IN) :: IND
		INTEGER, INTENT(IN) :: IB
		INTEGER, INTENT(IN) :: D123
		INTEGER, INTENT(IN) :: N
		INTEGER, INTENT(IN) :: MOP
		
		INTEGER, INTENT(IN) :: DIM_KP(SD),DIM_OP(SD)
		INTEGER, INTENT(IN) :: SPAN(NUMEL,SD) ! SPAN(IND,J) = I 
		
		DOUBLE PRECISION, INTENT(IN) :: BAK_KNOT0(NUMEL*D123*SD)!U0
		DOUBLE PRECISION, INTENT(IN) :: BAK_KNOT((DIM_KP(1)+1)+(DIM_KP(2)+1)+(DIM_KP(3)+1))!U
		
		DOUBLE PRECISION, INTENT(OUT) :: BSS_BASIS(NUMEL*D123, &
		   2*(DIM_OP(1)+1)+2*(DIM_OP(2)+1)+2*(DIM_OP(3)+1))
		   
		DOUBLE PRECISION, INTENT(OUT) :: BSS_A(NUMEL*D123,2*(MOP+1))
		DOUBLE PRECISION, INTENT(OUT) :: BSS_NDU(NUMEL*D123,(MOP+1)*(MOP+1))
		DOUBLE PRECISION, INTENT(OUT) :: LEFT(NUMEL*D123,MOP),RIGHT(NUMEL*D123,MOP)
		
		DOUBLE PRECISION :: D
		DOUBLE PRECISION :: SAVED,TEMP
		
		INTEGER :: OUTER
		INTEGER :: I,J,K,R
		INTEGER :: S1,S2
		INTEGER :: RK,PK
		INTEGER :: J1,J2
		
		INTEGER :: BSS_IND_1,BSS_IND_2
		INTEGER :: BAK_IND_1
		
		BSS_IND_1 = (IND-1)*D123 + IB
		BSS_IND_2 = 0
		BAK_IND_1 = 0
		
		DO OUTER = 1,SD
		  
		  BSS_NDU(BSS_IND_1,1) = 1.0D0
		  
		  DO J = 1,DIM_OP(OUTER)
		    
		    LEFT(BSS_IND_1,J) = BAK_KNOT0((OUTER-1)*NUMEL*D123+(IND-1)*D123+IB) - &
		       BAK_KNOT(BAK_IND_1+(SPAN(IND,OUTER)+1-J)+1)
		    RIGHT(BSS_IND_1,J) = BAK_KNOT(BAK_IND_1+(SPAN(IND,OUTER)+J)+1) - &
		       BAK_KNOT0((OUTER-1)*NUMEL*D123+(IND-1)*D123+IB)
		    SAVED = 0.0D0
		    
		    DO R = 0,J-1
		      
		      BSS_NDU(BSS_IND_1,J*(MOP+1)+R+1) = RIGHT(BSS_IND_1,R+1) + LEFT(BSS_IND_1,J-R)
		      TEMP = BSS_NDU(BSS_IND_1,R*(MOP+1)+J) / BSS_NDU(BSS_IND_1,J*(MOP+1)+R+1)
		      
		      BSS_NDU(BSS_IND_1,R*(MOP+1)+J+1) = SAVED + RIGHT(BSS_IND_1,R+1) * TEMP
		      SAVED = LEFT(BSS_IND_1,J-R) * TEMP
		      
		    END DO
		    
		    BSS_NDU(BSS_IND_1,J*(MOP+1)+J+1) = SAVED
		  
		  END DO	
		  
		  DO J = 0,DIM_OP(OUTER)
		    BSS_BASIS(BSS_IND_1,BSS_IND_2+J+1) = BSS_NDU(BSS_IND_1,J*(MOP+1)+DIM_OP(OUTER)+1)
		  END DO
		  BSS_IND_2 = BSS_IND_2 + (DIM_OP(OUTER)+1)
		  ! ----------------------------------------
		  
		  DO R = 0,DIM_OP(OUTER)
	  	    S1 = 0
	  	    S2 = (MOP+1)
	  	    BSS_A(BSS_IND_1,1) = 1.0D0
	  	    DO K = 1,N
	  	    
	  	      D = 0.0D0
	  	      RK = R-K
	  	      PK = DIM_OP(OUTER)-K
	  	      
	  	      IF (R .GE. K) THEN
	  	        BSS_A(BSS_IND_1,S2+1) = BSS_A(BSS_IND_1,S1+1) / &
	  	           BSS_NDU(BSS_IND_1,(PK+1)*(MOP+1)+RK+1)
		        D = BSS_A(BSS_IND_1,S2+1) * BSS_NDU(BSS_IND_1,RK*(MOP+1)+PK+1)
	  	      END IF
		        
		      IF (RK .GE. -1) THEN
		        J1 = 1
		      ELSE
		        J1 = -RK
		      END IF
		      
		      IF (R-1 .LE. PK) THEN
		        J2 = K-1
		      ELSE
		        J2 = DIM_OP(OUTER)-R
		      END IF
		      
		      DO J = J1,J2
		        BSS_A(BSS_IND_1,S2+J+1) = (BSS_A(BSS_IND_1,S1+J+1)-BSS_A(BSS_IND_1,S1+J)) / &
		           BSS_NDU(BSS_IND_1,(PK+1)*(MOP+1)+RK+J+1)
		        D = D + BSS_A(BSS_IND_1,S2+J+1) * BSS_NDU(BSS_IND_1,(RK+J)*(MOP+1)+PK+1)
		      END DO
		      
		      IF (R .LE. PK) THEN
		        BSS_A(BSS_IND_1,S2+K+1) = -BSS_A(BSS_IND_1,S1+K) / &
		           BSS_NDU(BSS_IND_1,(PK+1)*(MOP+1)+R+1)
		        D = D + BSS_A(BSS_IND_1,S2+K+1) * BSS_NDU(BSS_IND_1,R*(MOP+1)+PK+1)
		      END IF
		      
		      BSS_BASIS(BSS_IND_1,BSS_IND_2+R+1) = D
		      
		      J = S1
		      S1 = S2
		      S2 = J
		    
		    END DO
		  END DO
		  
		  R = DIM_OP(OUTER)
		  DO K = 1,N
		    DO J = 0,DIM_OP(OUTER)
		      BSS_BASIS(BSS_IND_1,BSS_IND_2+J+1) = &
		         BSS_BASIS(BSS_IND_1,BSS_IND_2+J+1) * DBLE(R)
		    END DO
		    R = R * (DIM_OP(OUTER)-K)
		  END DO
		  
		  BSS_IND_2 = BSS_IND_2 + (DIM_OP(OUTER)+1)
		  BAK_IND_1 = BAK_IND_1 + (DIM_KP(OUTER)+1)
		  
		END DO
		
	END SUBROUTINE DERS_BASIS_FUNS
	
	SUBROUTINE GRAD(SD,NEN,X,DPHI,FX)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD
		INTEGER, INTENT(IN) :: NEN
		DOUBLE PRECISION, INTENT(IN) :: DPHI(NEN,SD)
		DOUBLE PRECISION, INTENT(IN) :: X(NEN,SD)
		DOUBLE PRECISION, INTENT(OUT) :: FX(SD,SD)
		
		INTEGER :: J1,J2,K
		
		FX(:,:) = 0.0D0
		
		DO J1=1,SD
		  DO J2=1,SD
		    DO K = 1,NEN
		      FX(J1,J2) = FX(J1,J2) + X(K,J1)*DPHI(K,J2)
		    END DO
		  END DO
		END DO
		
	END SUBROUTINE GRAD
	
	SUBROUTINE INV(SD,F,INVF,DETF)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD
		DOUBLE PRECISION, INTENT(IN) :: F(SD,SD)
		DOUBLE PRECISION, INTENT(OUT) :: INVF(SD,SD)
		DOUBLE PRECISION, INTENT(OUT) :: DETF
		
		DOUBLE PRECISION :: ADJF(SD,SD)
		
		INTEGER :: J1,J2
		
		CALL ADJ(SD,F,ADJF)
		DETF = DET(SD,F)
		DO J1 = 1,SD
		  DO J2 = 1,SD
		    INVF(J1,J2) = ADJF(J2,J1) / DETF
		  END DO
		END DO
		
	END SUBROUTINE INV
	
	SUBROUTINE ADJ(SD,F,ADJF)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD
		DOUBLE PRECISION, INTENT(IN) :: F(SD,SD)
		DOUBLE PRECISION, INTENT(OUT) :: ADJF(SD,SD)
		
		IF (SD .EQ. 1) THEN
		  ADJF(1,1) = 1.0D0
		ELSE IF (SD .EQ. 2) THEN
		  ADJF(1,1) =  F(2,2)
		  ADJF(1,2) = -F(2,1)
		  ADJF(2,1) = -F(1,2)
		  ADJF(2,2) =  F(1,1)
		ELSE IF (SD .EQ. 3) THEN
		  ADJF(1,1) = F(2,2)*F(3,3)-F(3,2)*F(2,3)
		  ADJF(1,2) = F(3,1)*F(2,3)-F(2,1)*F(3,3)
		  ADJF(1,3) = F(2,1)*F(3,2)-F(3,1)*F(2,2)
		  ADJF(2,1) = F(1,3)*F(3,2)-F(1,2)*F(3,3)
		  ADJF(2,2) = F(1,1)*F(3,3)-F(1,3)*F(3,1)
		  ADJF(2,3) = F(1,2)*F(3,1)-F(1,1)*F(3,2)
		  ADJF(3,1) = F(1,2)*F(2,3)-F(1,3)*F(2,2)
		  ADJF(3,2) = F(1,3)*F(2,1)-F(2,3)*F(1,1)
		  ADJF(3,3) = F(1,1)*F(2,2)-F(2,1)*F(1,2)
		END IF
		
	END SUBROUTINE ADJ
	
	DOUBLE PRECISION FUNCTION DET(SD,F)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: SD
		DOUBLE PRECISION, INTENT(IN) :: F(SD,SD)	
		
		IF (SD .EQ. 1) THEN
		  DET = F(1,1)
		ELSE IF (SD .EQ. 2) THEN
		  DET = F(1,1)*F(2,2) - F(1,2)*F(2,1)
		ELSE IF (SD .EQ. 3) THEN
		  DET = F(1,1)*( F(2,2)*F(3,3)-F(3,2)*F(2,3) ) &
		      + F(1,2)*( F(3,1)*F(2,3)-F(2,1)*F(3,3) ) &
		      + F(1,3)*( F(2,1)*F(3,2)-F(3,1)*F(2,2) )
		END IF
		
	END FUNCTION DET
	
	SUBROUTINE MAT_TRANS(A,N,M,AT)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: N, M
		DOUBLE PRECISION, INTENT(IN) :: A(N,M)
		DOUBLE PRECISION, INTENT(OUT) :: AT(M,N)
		
		INTEGER :: J1,J2
		
		DO J1 = 1,N
		  DO J2 = 1,M
		    AT(J2,J1) = A(J1,J2)
		  END DO
		END DO
		
	END SUBROUTINE MAT_TRANS
	
	SUBROUTINE MATVEC_PRO(A,N,M,B,C)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: N, M
		DOUBLE PRECISION, INTENT(IN) :: A(N,M), B(M)
		DOUBLE PRECISION, INTENT(OUT) :: C(N)
		
		INTEGER :: J,K
		
		C(:) = 0.0D0
		DO J = 1, N
		  DO K = 1, M
		    C(J) = C(J) + A(J,K)*B(K)
		  END DO
		END DO
		
	END SUBROUTINE MATVEC_PRO
	

END MODULE M_CPU
