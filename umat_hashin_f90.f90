!     ------------------------------------------------------------------
!     |     DATE           PROGRAMMER           DESCRIPTION OF CHANGE  |  
!     |     ====           ==========            ===================== |
!     |  29/09/2020      MÁRIO RUI ARRUDA                              |
!     |  IST LISBON                                                    |
!     ------------------------------------------------------------------

!     ----------UMAT ONLY FOR PLANE ELEMENST WITH HASHIN DAMAGE---------

!     #COPYRIGHT 2020 BY MÁRIO RUI TIAGO ARRUDA
!     ALL RIGHTS RESERVED. NO PART OF THIS SUBROUTINE MAY BE REPRODUCED OR 
!     USED IN ANY MANNER WITHOUT WRITTEN PERMISSION OF THE COPYRIGHT OWNER.

! When using this UMAT in papers, always cite these works https://doi.org/10.1016/j.compositesb.2020.107818 and 
! https://doi.org/10.1016/j.compstruct.2020.112453 and always provide credit to the original authors

!     ------------------------------------------------------------------
!     ----------ABAQUS INPUT VARIABLES IN UMAT SUBROUTINE---------------
!     ------------------------------------------------------------------

subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, &
           ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, & 
           nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, &
           npt, layer, kspt, kstep, kinc)

IMPLICIT NONE
 

!     ------------------------------------------------------------------
!     -------------------ABAQUS DIMENSION VARIABLES---------------------
!     ------------------------------------------------------------------
 
!     -----------------ABAQUS CHARACTER VARIABLES----------------------
character(len=80) :: cmname 

!     -------------------ABAQUS INTEGER VARIABLES-----------------------  
integer :: ntens,ndi,nshr,nstatv,nprops,noel,npt,kspt,kstep,kinc,nprecd,layer

!     -------------------ABAQUS REAL VARIABLES--------------------------
real(kind=8) :: celent,sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt 

!     -------------------ABAQUS ARRAY VARIABLES-------------------------
real(kind=8) :: stress(ntens),statev(nstatv),ddsdde(ntens,ntens),&
                ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),&
                predef(1),dpred(1),props(nprops),coords(3),drot(3,3),&
                dfgrd0(3,3),dfgrd1(3,3),time(2)

					 
!     ------------------------------------------------------------------    
!     -----------------DECLARATION OF VARIABLES-------------------------
!     ------------------------------------------------------------------

integer :: i,j,k,l,m,n,kiter,ktotal

real(kind=8) :: E1,E2,G,v12,v21,Xt,Xc,Yt,Yc,St,Sl,Faft,Fafc,Famt,Famc,&
								afto,afco,amto,amco,Gft,Gfc,Gmt,Gmc,eta,dmax,dfto,dfco,dmto,dmco,dso,&
								dvft,dvfc,dvmt,dvmc,dvfto,dvfco,dvmto,dvmco,dvso,seqft0,seqfc0,seqmt0,seqmc0,&
								df,dvf,dm,dvm,Dcoef,catLc,dft,dfc,dmt,dmc,ds,dvs,Fafto,Fafco,Famto,Famco,&
								ueqft,ueqft0,ueqftu,ueqfc,ueqfc0,ueqfcu,ueqmt,ueqmt0,ueqmtu,ueqmc,ueqmc0,ueqmcu,&
                ueqftp,ueqfcp,ueqmtp,ueqmcp,pft,pfc,pmt,pmc,alpha

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4

!     -------------INITIATIONS OF ARRAYS--------------------------------
real(kind=8) :: eij(ntens),eijbp(ntens),eijbn(ntens),eijo(ntens),eijbpo(ntens),eijbno(ntens),&
                sij(ntens),sijbp(ntens),sijbn(ntens),sijo(ntens),sijbpo(ntens),sijbno(ntens),&
                sije(ntens),sijebp(ntens),sijebn(ntens),sijeo(ntens),sijebpo(ntens),sijebno(ntens),&
                ddsddee(ntens,ntens),eij_e(ntens)
                
!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3-------------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0, FOUR=4.d0,&
                           HALF=0.5d0,TOL0=1.d-8, TOLseq=1.d-2			
						
!     --------------ELASTIC AND MECHANICAL PROPERTIES-------------------

E1=props(1)     ! LONGITUDINAL ELASTIC MODULUS
E2=props(2)     ! TRANSVERSAL ELASTIC MODULUS
G=props(3)      ! SHEAR ELASTIC MODULUS
v12=props(4)    ! POISSON COEFICIENT
v21=v12*E2/E1   ! SYMMETRIC POISSON COEFICIENT

Xt=props(5)     ! LONGITUDINAL TENSILE RESISTENCE
Xc=props(6)     ! LONGITUDINAL COMPRESSIVE RESISTENCE
Yt=props(7)     ! TRANSVERSAL TENSILE RESISTENCE
Yc=props(8)     ! TRANSVERSAL COMPRESSIVE RESISTENCE
Sl=props(9)     ! LONGITUDINAL SHEAR RESISTENCE
St=props(10)    ! TRANSVERSAL SHEAR RESISTENCE
alpha=props(11)

Gft=props(12)   ! FRACTURE ENERGY FOR FIBRE TENSION
Gfc=props(13)   ! FRACTURE ENERGY FOR FIBRE COMPRESSION
Gmt=props(14)   ! FRACTURE ENERGY FOR MATRIX TENSION
Gmc=props(15)   ! FRACTURE ENERGY FOR MATRIX COMPRESSION

eta=props(16)   ! VISCOUS REGULARIZATION COEFICIENT
dmax=props(17)  ! MAXIMUM ALLOWED DAMAGE FOR CONVERGENCY

pft=props(18)   ! RESIDUAL LONGITUDINAL TENSILE RESISTENCE
pfc=props(19)   ! RESIDUAL LONGITUDINAL COMPRESSIVE RESISTENCE
pmt=props(20)   ! RESIDUAL TRANSVERSAL TENSILE RESISTENCE
pmc=props(21)   ! RESIDUAL TRANSVERSAL COMPRESSIVE RESISTENCE

!     ---------------STATE FIELD VARIABLES FOR ABAQUS-------------------
if (kinc == 1) then 
	do i=1, nstatv
		statev(i)=ZERO  ! VARIABLES INITIATION (if not=0, depending on compiler LINUX/WINDOWS)
  end do
end if

dfto=statev(1)      ! FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dfco=statev(2)      ! FIBRE COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dmto=statev(3)      ! MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dmco=statev(4)      ! MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dso=statev(5)       ! SHEAR DAMAGE FROM PREVIOUS INCREMENT

Fafto=statev(6)     ! FIBRE TENSION CRITERIA FROM PREVIOUS INCREMENT
Fafco=statev(7)     ! FIBRE COMPRESSION CRITERIA FROM PREVIOUS INCREMENT
Famto=statev(8)     ! MATRIX TENSION CRITERIA FROM PREVIOUS INCREMENT
Famco=statev(9)     ! MATRIX COMPRESSION CRITERIA FROM PREVIOUS INCREMENT

dvfto=statev(10)    ! VISCOUS FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dvfco=statev(11)    ! VISCOUS FIBRE COMRESSION DAMAGE FROM PREVIOUS INCREMENT
dvmto=statev(12)    ! VISCOUS MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dvmco=statev(13)    ! VISCOUS MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dvso=statev(14)     ! VISCOUS SHEAR DAMAGE FROM PREVIOUS INCREMENT

seqft0=statev(15)   ! EQUIVALENTE FIBRE TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqfc0=statev(16)   ! EQUIVALENTE FIBRE COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmt0=statev(17)   ! EQUIVALENTE MATRIX TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmc0=statev(18)   ! EQUIVALENTE MATRIX COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

ueqft0=statev(19)   ! EQUIVALENTE FIBRE TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqfc0=statev(20)   ! EQUIVALENTE FIBRE COMPRESSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmt0=statev(21)   ! EQUIVALENTE MATRIX TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmc0=statev(22)   ! EQUIVALENTE MATRIX COMPRESSON DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

sijeo(1)=statev(23) ! EFFECTIVE S11 STRESS AT THE BEGGING OF THE INCREMENT
sijeo(2)=statev(24) ! EFFECTIVE S22 STRESS AT THE BEGGING OF THE INCREMENT
sijeo(3)=statev(25) ! EFFECTIVE S12 STRESS AT THE BEGGING OF THE INCREMENT


!     ------------------------------------------------------------------
!     ---------FIELD PREDICTOR FOR VARIABLES AND STIFFNESS -------------
!     ------------------------------------------------------------------

!     --------------FINAL STRAIN INCREMENT------------------------------
do i=1, ntens
  eij(i)=stran(i)+dstran(i) ! TOTAL STRAIN IN THE BEGGIND OF THE INCREMENT
  eijo(i)=stran(i)          ! STRAIN FROM PREVIOUS INCREMENT
  sijo(i)=stress(i)         ! STRESS FROM PREVIOUS INCREMENT
end do

!     ----STABILITY CONDITION DURING EXPLICIT MATERIAL INTEGRATION------
do i=1,ntens
	if (abs(dstran(i)) > (Xt/E1)) then
		pnewdt=0.5
	end if
end do

!     -------VERIFICATION OF FIBRE AND MATRIX TENSION/COMPRESSION-------
if (sijeo(1) >= ZERO) then ! FIBRE IN TENSION VS COMPRESSION
	df=dfto
	dvf=dvfto
else
	df=dfco
	dvf=dvfco
end if   
if (sijeo(2) >= ZERO) then ! MATRIX IN TENSION VS COMPRESSION
	dm=dmto
	dvm=dvmto
else
	dm=dmco
	dvm=dvmco
end if

ds=dso

Faft=Fafto ! INITIATION FIBRE TENSION FAILURE
Fafc=Fafco ! INITIATION FIBRE COMPRESSION FAILURE
Famt=Famto ! INITIATION MATRIX TENSION FAILURE
Famc=Famco ! INITIATION MATRIX COMPRESSION FAILURE

!     -----SECANT DAMAGED MATRIX AT THE BEGINNING OF THE INCREMENT------
call calc_stiff_assembler(ntens,ddsdde,Dcoef,E1,E2,G,df,dm,ds,v12,v21)

!     -----SECANT EFFECTIVE MATRIX AT THE BEGINNING OF THE INCREMENT----
call calc_efect_stiff_assembler(ntens,ddsddee,Dcoef,E1,E2,G,df,dm,ds,v12,v21)

!    -----------MULTIPLICATION OF SECANT MATRIX WITH STRAIN-------------
call calc_s_ke_multiply(ndi,ntens,sij,eij,ddsdde)
call calc_s_ke_multiply(ndi,ntens,sije,eij,ddsddee)

!     -------------ZERO STRESS IN CASE OF PURE SHEAR STRESS-------------
if (abs(sije(1)) <= (ZERO+TOL0)) then
	sije(1)=ZERO
endif
if (abs(sije(2)) <= (ZERO+TOL0)) then
	sije(2)=ZERO
endif 


!     ------------------------------------------------------------------
!     ------------------CREATING MACAULY BRACKETS-----------------------
!     ------------------------------------------------------------------
call calc_macaulay_bracket(ndi,ntens,eij,eijbp,eijbn)        ! BRACKETS FOR STRAIN
call calc_macaulay_bracket(ndi,ntens,sij,sijbp,sijbn)        ! BRACKETS FOR STRESS
call calc_macaulay_bracket(ndi,ntens,sije,sijebp,sijebn)     ! BRACKETS FOR EFECTIVE STRESS
call calc_macaulay_bracket(ndi,ntens,eijo,eijbpo,eijbno)     ! BRACKETS FOR PREVIOUS STRAIN
call calc_macaulay_bracket(ndi,ntens,sijo,sijbpo,sijbno)     ! BRACKETS FOR PREVIOUS STRESS
call calc_macaulay_bracket(ndi,ntens,sijeo,sijebpo,sijebno)  ! BRACKETS FOR PREVIOUS EFECTIVE STRESS


!     ------------------------------------------------------------------
!     --------------FAILURE CRITERIA VERIFICATION-----------------------
!     ------------------------------------------------------------------

!     -------------CARACTERISTIC LENGTH FROM ABAQUS---------------------
catLc=celent

!     --------VERIFICATION OF FIBRE TENSION/COMPRESSION FAILURE---------

if (sije(1) >= ZERO .and. Faft < ONE) then
  Faft=((sije(1)/Xt)**2+alpha*(sije(3)/Sl)**2)
  if (Faft >= ONE) then
    call calc_activation(catLc,Fafto,alpha,eijbpo(1),eijo(3),sijbpo(1),sijo(3),ueqft0,seqft0)
    statev(19)=ueqft0
    statev(15)=seqft0
    Faft=ONE
    statev(6)=Faft
  end if
else if (sije(1) < ZERO .and. Fafc < ONE) then
  Fafc=(sije(1)/Xc)**2
	if (Fafc >= ONE) then
    call calc_activation(catLc,Fafco,ZERO,eijbno(1),eijo(3),sijbno(1),sijo(3),ueqfc0,seqfc0)
    statev(20)=ueqfc0
    statev(16)=seqfc0
    Fafc=ONE
    statev(7)=Fafc
	end if
end if

!     --------VERIFICATION OF MATRIX TENSION/COMPRESSION FAILURE---------
if (sije(2) >= ZERO .and. Famt < ONE) then
  Famt=((sije(2)/Yt)**2+(sije(3)/Sl)**2)
  if (Famt >= ONE) then
    call calc_activation(catLc,Famto,ONE,eijbpo(2),eijo(3),sijbpo(2),sijo(3),ueqmt0,seqmt0)
    statev(21)=ueqmt0
    statev(17)=seqmt0
    Famt=ONE
    statev(8)=Famt
  end if
else if (sije(2) < ZERO .and. Famc < ONE) then
  Famc=((HALF*sije(2)/St)**2+((HALF*Yc/St)**2-ONE)*sije(2)/Yc+(sije(3)/Sl)**2)
	if (Famc >= ONE) then
    call calc_activation_extra(catLc,Yc,St,Sl,sijeo(2),sijeo(3),eijbno(2),eijo(3),sijbno(2),sijo(3),ueqmc0,seqmc0)
    statev(22)=ueqmc0
    statev(18)=seqmc0
    Famc=ONE
    statev(9)=Famc
	end if
end if

statev(6)=Faft
statev(7)=Fafc
statev(8)=Famt
statev(9)=max(ZERO,Famc)


!     ------------------------------------------------------------------
!     -----------------HASHIN DAMAGE EVOLUTION--------------------------
!     ------------------------------------------------------------------

dft=dfto;dfc=dfco;dmt=dmto;dmc=dmco;
dvft=dvfto;dvfc=dvfco;dvmt=dvmto;dvmc=dvmco;

!     -------------FIBER TENSION DAMAGE EVOLUTION-----------------------
if (Faft >= ONE) then
	ueqft=sqrt(eijbp(1)**2+alpha*eij(3)**2)*catLc
  call calc_damage_lin(ueqft,ueqft0,ueqftu,eta,dft,dfto,dvft,dvfto,dmax,dtime,pft,ueqftp,seqft0,Gft)
  !call calc_damage_exp(ueqft,ueqft0,ueqftu,eta,dft,dfto,dvft,dvfto,dmax,dtime,pft,ueqftp,seqft0,Gft)
  statev(1)=dft
	statev(10)=dvft
end if

!     -------------FIBER COMPRESSION DAMAGE EVOLUTION-------------------
if (Fafc >= ONE) then
	ueqfc=sqrt(eijbn(1)**2)*catLc
  call calc_damage_lin(ueqfc,ueqfc0,ueqfcu,eta,dfc,dfco,dvfc,dvfco,dmax,dtime,pfc,ueqfcp,seqfc0,Gfc)
  !call calc_damage_exp(ueqfc,ueqfc0,ueqfcu,eta,dfc,dfco,dvfc,dvfco,dmax,dtime,pfc,ueqfcp,seqfc0,Gfc)
  statev(2)=dfc
  statev(11)=dvfc
end if

!     -------------MATRIX TENSION DAMAGE EVOLUTION----------------------
if (Famt >= ONE) then
	ueqmt=sqrt(eijbp(2)**2+eij(3)**2)*catLc
  call calc_damage_lin(ueqmt,ueqmt0,ueqmtu,eta,dmt,dmto,dvmt,dvmto,dmax,dtime,pmt,ueqmtp,seqmt0,Gmt)
  !call calc_damage_exp(ueqmt,ueqmt0,ueqmtu,eta,dmt,dmto,dvmt,dvmto,dmax,dtime,pmt,ueqmtp,seqmt0,Gmt)
  statev(3)=dmt
	statev(12)=dvmt
end if

!     -----------MATRIX COMPRESSION DAMAGE EVOLUTION--------------------
if (Famc >= ONE) then
	ueqmc=sqrt(eijbn(2)**2+eij(3)**2)*catLc
  call calc_damage_lin(ueqmc,ueqmc0,ueqmcu,eta,dmc,dmco,dvmc,dvmco,dmax,dtime,pmc,ueqmcp,seqmc0,Gmc)
  !call calc_damage_exp(ueqmc,ueqmc0,ueqmcu,eta,dmc,dmco,dvmc,dvmco,dmax,dtime,pmc,ueqmcp,seqmc0,Gmc)
  statev(4)=dmc
	statev(13)=dvmc
end if

!     --------------------SHEAR DAMAGE EVOLUTION------------------------
ds=ONE-(ONE-dft)*(ONE-dfc)*(ONE-dmt)*(ONE-dmc)
dvs=ONE-(ONE-dvft)*(ONE-dvfc)*(ONE-dvmt)*(ONE-dvmc)
statev(5)=ds
statev(14)=dvs


!     ------------------------------------------------------------------
!     -----------------FINAL DAMAGE BEHAVIOUR---------------------------
!     ------------------------------------------------------------------

!     -----------------------FIBRE DAMAGE-------------------------------
if (sije(1) >= ZERO) then
	df=dft
	dvf=dvft
else
	df=dfc
	dvf=dvfc
endif
		
!     -----------------------MATRIX DAMAGE------------------------------  
if (sije(2) >= ZERO) then
	dm=dmt
	dvm=dvmt
else
	dm=dmc
	dvm=dvmc
endif


!     ------------------------------------------------------------------
!     -------FIELD CORRECTOR FOR STRESS AND SECANT STIFFNESS -----------
!     ------------------------------------------------------------------

!     ---------VISCOUS SECANT DAMAGED MATRIX IN THE ITERATION-----------
call calc_stiff_assembler(ntens,ddsdde,Dcoef,E1,E2,G,dvf,dvm,dvs,v12,v21)

!     --------VISCOUS EFFECTIVE SECANT MATRIX IN THE ITERATION----------
call calc_efect_stiff_assembler(ntens,ddsddee,Dcoef,E1,E2,G,dvf,dvm,dvs,v12,v21)

!    ---------------------VISCOUS STRESS CORRECTOR----------------------
call calc_s_ke_multiply(ndi,ntens,stress,eij,ddsdde)
call calc_s_ke_multiply(ndi,ntens,sije,eij,ddsddee)

statev(23)=sije(1)
statev(24)=sije(2)
statev(25)=sije(3)

!    -----------------UPDATED DEFORMATION ENERGY------------------------
call calc_half_s_e_multiply(ndi,ntens,sse,stress,sijo,dstran)

return
end


!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||USER SUBROUTINES||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


!     __________________________________________________________________
!     ___________SUBROUTINE FOR MACAULY BRACKET OPERATION_______________
!     __________________________________________________________________

subroutine calc_macaulay_bracket(ndi,ntens,mat,mat_bp,mat_bn)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: mat(ntens)
real(kind=8), intent(inout) :: mat_bp(ntens),mat_bn(ntens)

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2---------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
	mat_bp(i)=HALF*(mat(i)+abs(mat(i)))  ! BRACKETS FOR DAMAGE
	mat_bn(i)=HALF*(-mat(i)+abs(mat(i))) ! BRACKETS FOR DAMAGE
end do

return
end subroutine


!     __________________________________________________________________
!     _______SUBROUTINE FOR STRAIN STIFFNESS MULTIPLICATION_____________
!     __________________________________________________________________

subroutine calc_s_ke_multiply(ndi,ntens,mat_si,mat_ej,mat_Cij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: mat_Cij(ntens,ntens), mat_ej(ntens)
real(kind=8), intent(inout) :: mat_si(ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
  mat_si(i)=ZERO
  do j=1,ndi ! NORMAL STRAIN
	  mat_si(i)=mat_si(i)+mat_Cij(i,j)*mat_ej(j)
  enddo
enddo
do i=ndi+1,ntens ! SHEAR STRAIN
  mat_si(i)=mat_Cij(i,i)*mat_ej(i) 
enddo

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FOR DEFORMATION ENERGY__________________
!     __________________________________________________________________

subroutine calc_half_s_e_multiply(ndi,ntens,energy,sij,sij_old,deij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: sij(ntens), sij_old(ntens), deij(ntens)
real(kind=8), intent(inout) :: energy

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
  energy=energy+HALF*(sij_old(i)+sij(i))*deij(i)
end do

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE STIFFNESS ASSEMBLER____________________
!     __________________________________________________________________

subroutine calc_stiff_assembler(ntens,Cij,Dvar,E1,E2,G,df,dm,ds,v12,v21)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens

!     ----------------------REAL MATRIX VARIABLES----------------------- 
real(kind=8), intent(in)  :: E1,E2,G,df,dm,ds,v12,v21
real(kind=8), intent(inout) :: Cij(ntens,ntens), Dvar

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21    
Cij(1,1)=(ONE-df)*E1*ONE/Dvar
Cij(2,2)=(ONE-dm)*E2*ONE/Dvar
Cij(3,3)=(ONE-ds)*G    
Cij(1,2)=(ONE-df)*(ONE-dm)*v21*E1*ONE/Dvar
Cij(2,1)=(ONE-df)*(ONE-dm)*v12*E2*ONE/Dvar

return
end subroutine


!     __________________________________________________________________
!     ____________SUBROUTINE EFECTIVE STIFFNESS ASSEMBLER_______________
!     __________________________________________________________________

subroutine calc_efect_stiff_assembler(ntens,Cij,Dvar,E1,E2,G,df,dm,ds,v12,v21)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: E1,E2,G,df,dm,ds,v12,v21
real(kind=8), intent(inout) :: Cij(ntens,ntens), Dvar

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21    
Cij(1,1)=E1*ONE/Dvar
Cij(2,2)=E2*ONE/Dvar
Cij(3,3)=G    
Cij(1,2)=(ONE-dm)*v21*E1*ONE/Dvar
Cij(2,1)=(ONE-df)*v12*E2*ONE/Dvar

return
end subroutine


!     __________________________________________________________________
!     ______________SUBROUTINE FOR ACTIVATION FUNCTION__________________
!     __________________________________________________________________

subroutine calc_activation(catLc,Fa,alpha,eii,eij,sii,sij,ueq0,seq0)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8) :: root
real(kind=8), intent(in)  :: catLc,Fa,alpha,eii,eij,sii,sij
real(kind=8), intent(inout) :: ueq0,seq0

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

root=ONE/sqrt(Fa)
ueq0=root*sqrt(eii**2+alpha*eij**2)*catLc
seq0=root**2*(sii*eii+alpha*sij*eij)/(ueq0/catLc)

return
end subroutine


!     __________________________________________________________________
!     ___________SUBROUTINE FOR EXTRA ACTIVATION FUNCTION_______________
!     __________________________________________________________________

subroutine calc_activation_extra(catLc,Yc,St,Sl,siie,sije,eii,eij,sii,sij,ueq0,seq0)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8) :: a,b,c,root
real(kind=8), intent(in)  :: catLc,siie,sije,eii,eij,sii,sij,Yc,St,Sl
real(kind=8), intent(inout) :: ueq0,seq0

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

a=(HALF*siie/St)**2+(sije/Sl)**2
b=((HALF*Yc/St)**2-ONE)*siie/Yc
c=-ONE
root=abs((-b+sqrt(b**2-FOUR*a*c))/(TWO*a))
ueq0=root*sqrt(eii**2+eij**2)*catLc
seq0=root**2*(sii*eii+sij*eij)/(ueq0/catLc)

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE FOR DAMAGE EVOLUTION___________________
!     __________________________________________________________________

subroutine calc_damage_lin(ueq,ueq0,uequ,eta,d,dold,dv,dvold,dmax,dtime,pres,ueqp,seq0,G)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8), intent(in)  :: ueq,ueq0,eta,dmax,dvold,dold,dtime,pres,seq0,G
real(kind=8), intent(inout) :: d,dv,ueqp,uequ

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

uequ=ueq0+TWO*G/seq0
ueqp=uequ-pres*(uequ-ueq0)
if (ueq > ueq0 .and. ueq<ueqp) then
  d=uequ*(ueq-ueq0)/(ueq*(uequ-ueq0))
	dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
else if (ueq>=ueqp) then
  d=ONE-pres*ueq0/ueq
  dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
end if

d=min(d,dmax)
d=max(d,dold)
dv=min(dv,dmax)
dv=max(dv,dvold)

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE FOR DAMAGE EVOLUTION___________________
!     __________________________________________________________________

subroutine calc_damage_exp(ueq,ueq0,uequ,eta,d,dold,dv,dvold,dmax,dtime,pres,ueqp,seq0,G)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8), intent(in)  :: ueq,ueq0,eta,dmax,dvold,dold,dtime,pres,seq0,G
real(kind=8), intent(inout) :: d,dv,ueqp,uequ

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (pres==ZERO) then
  if (ueq > ueq0) then
    d=ONE-ueq0/ueq*exp(-seq0/G*(ueq-ueq0))
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  end if
else
  ueqp=-G/seq0*log(pres)+ueq0
  if (ueq > ueq0 .and. ueq<ueqp) then
    d=ONE-ueq0/ueq*exp(-seq0/G*(ueq-ueq0))
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  else if (ueq>=ueqp) then
    d=ONE-pres*ueq0/ueq
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  end if
end if

d=min(d,dmax)
d=max(d,dold)
dv=min(dv,dmax)
dv=max(dv,dvold)

return
end subroutine
