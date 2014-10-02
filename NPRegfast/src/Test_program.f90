	

program hola
implicit none
!integer,parameter::kbin=100,ncmax=5,nboot=10
integer i,p,nf,fact(100),iend,n,ikernel,iopt,model,kernel,nh
integer,allocatable::F(:)
integer,allocatable::nc(:)
double precision,allocatable:: X(:),Y(:),w(:),h(:),C2(:,:),Xb(:),&
pb(:,:,:),li(:,:,:),ls(:,:,:),Dif(:,:,:,:),Difi(:,:,:,:),Difs(:,:,:,:),&
difc(:,:,:),difci(:,:,:),difcs(:,:,:),pboot(:,:,:,:),&
c(:,:),cs(:,:),ci(:,:),pcmin(:),pcmax(:),Cboot(:,:,:),a(:),ainf(:),asup(:),b(:),binf(:),bsup(:)
double precision T,pvalor
integer kbin,ncmax,nboot
character*30 fichero





kbin=100
ncmax=5
nboot=100




!write (*,*) 'nombre del archivo de datos'
!read (*,*) fichero
!fichero="alo.txt"
!fichero="decus_phill_1.txt"
!fichero=".prn"
!fichero="pediatria3.txt"
fichero="Datos_a?o.prn"
n=0
open (1,file=fichero)
read (1,*,iostat=iend)
1   continue
    read (1,*,iostat=iend)
if (iend.eq.0) then
n=n+1
goto 1
end if
close(1)



allocate (F(n),X(n),Y(n),W(n),h(100))

W=1
open (1,file=fichero)
read (1,*)
do i=1,n
read(1,*) F(i),X(i),Y(i)
end do
close (1)





call Factores(F,n,fact,nf)

allocate (C2(ncmax,nf),nc(nf),xb(kbin),pb(kbin,3,nf),li(kbin,3,nf),&
ls(kbin,3,nf),Dif(kbin,3,nf,nf),Difi(kbin,3,nf,nf),Difs(kbin,3,nf,nf),&
difc(3,nf,nf),difci(3,nf,nf),difcs(3,nf,nf),pboot(kbin,3,nf,nboot),c(3,nf),cs(3,nf),ci(3,nf),pcmin(nf),pcmax(nf),Cboot(3,nf,nboot),&
a(nf),ainf(nf),asup(nf),b(nf),binf(nf),bsup(nf))
p=1
C2=-1
iopt=1
ikernel=1
h=-1
nc=-1

if(ikernel.eq.1) p=1

!h=98 !80
pcmax=200
pcmin=0



!model=1 noparametrico
!model=2 alometrico
!model=0 contraste

model=2
!loc=0
h=-1

kernel=1

nh=30
call frfast(F,X,Y,W,n,h,C2,nc,ncmax,p,kbin,fact,nf,ikernel,iopt,nboot,&
xb,pb,li,ls,dif,difi,difs,model,pvalor,c,cs,ci,difc,difcs,difci,T,pboot,pcmin,pcmax,cboot,kernel,nh,a,ainf,asup,b,binf,bsup)



write (*,*) (nc(i),i=1,nf)
write (*,*) (h(i),i=1,nf)

end
