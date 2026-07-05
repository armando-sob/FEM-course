!=============================================================
! Portico 2D - Fortran em formato livre
! GDL por no: ux, uy, rz
! Elemento: portico plano 2 nos, com rigidez elastica e P-Delta.
!=============================================================
program portico2d_pdelta
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 307)
   integer, parameter :: maxn=200, maxe=400, maxm=50, maxbc=600, maxl=600
   integer, parameter :: maxdof = 3*maxn

   integer :: nn, ne, nm, nbc, nload, method, ipdelta, maxiter
   integer :: ndof, nfree, iter, converged
   integer :: i, j, e, id, m, node, n1, n2, matid
   integer :: bcnode(maxbc), bcdof(maxbc), ldnode(maxl), lddof(maxl)
   integer :: e1(maxe), e2(maxe), emat(maxe), prescribed(maxdof)
   integer :: free_gdls(maxdof), g(6)
   real(dp) :: bcval(maxbc), ldval(maxl), x(maxn), y(maxn)
   real(dp) :: ea(maxm), aa(maxm), ei(maxm)
   real(dp) :: kmat(maxdof,maxdof), force(maxdof), disp(maxdof)
   real(dp) :: reaction(maxdof), disp_old(maxdof)
   real(dp) :: kwork(maxdof,maxdof), fwork(maxdof), ured(maxdof)
   real(dp) :: nax(maxe), nax_new(maxe)
   real(dp) :: kl(6,6), trans(6,6), ue(6), qel(6)
   real(dp) :: length, c, s, emod, area, iner, tol
   real(dp) :: err, den, diff, absu
   character(len=*), parameter :: infile='portico2d_pdelta.dat'
   character(len=*), parameter :: outfile='portico2d_pdelta.out'

   open(10, file=infile, status='old', action='read')
   open(20, file=outfile, status='replace', action='write')

   read(10, *) nn, ne, nm, nbc, nload
   read(10, *) method
   read(10, *) ipdelta, maxiter, tol
   ndof = 3*nn

   x=0.0_dp; y=0.0_dp; ea=0.0_dp; aa=0.0_dp; ei=0.0_dp
   kmat=0.0_dp; force=0.0_dp; disp=0.0_dp; disp_old=0.0_dp
   reaction=0.0_dp; prescribed=0; nax=0.0_dp; nax_new=0.0_dp

   do i=1,nn
      read(10, *) node, x(node), y(node)
   end do

   do i=1,nm
      read(10, *) m, ea(m), aa(m), ei(m)
   end do

   do e=1,ne
      read(10, *) id, e1(id), e2(id), emat(id)
   end do

   do i=1,nbc
      read(10, *) bcnode(i), bcdof(i), bcval(i)
      prescribed(gdl(bcnode(i), bcdof(i))) = 1
   end do

   do i=1,nload
      read(10, *) ldnode(i), lddof(i), ldval(i)
      force(gdl(ldnode(i), lddof(i))) = force(gdl(ldnode(i), lddof(i))) + ldval(i)
   end do
   close(10)

   converged = 0
   if (ipdelta == 0) then
      call build_stiffness(.false.)
      call solve_with_bc()
      call compute_axials(disp, nax)
      iter = 1
      converged = 1
   else
      if (maxiter <= 0) maxiter = 30
      if (tol <= 0.0_dp) tol = 1.0e-8_dp
      do iter=1,maxiter
         call build_stiffness(.true.)
         call solve_with_bc()
         call compute_axials(disp, nax_new)
         err = 0.0_dp
         den = 1.0_dp
         do i=1,ndof
            diff = abs(disp(i)-disp_old(i))
            absu = abs(disp(i))
            if (diff > err) err = diff
            if (absu > den) den = absu
         end do
         err = err/den
         disp_old = disp
         nax = nax_new
         if (iter > 1 .and. err <= tol) then
            converged = 1
            exit
         end if
      end do
   end if

   call build_stiffness(ipdelta == 1)
   reaction(1:ndof) = matmul(kmat(1:ndof,1:ndof), disp(1:ndof)) - force(1:ndof)

   write(20, *) '========================================'
   write(20, *) 'PORTICO 2D - RESULTADOS'
   write(20, *) 'Arquivo de entrada: ', infile
   write(20, *) 'NN, NE = ', nn, ne
   if (method == 1) then
      write(20, *) 'Metodo: PENALIDADE'
   else
      write(20, *) 'Metodo: ELIMINACAO DE LINHAS/COLUNAS'
   end if
   if (ipdelta == 1) then
      write(20, *) 'Analise: P-DELTA'
      write(20, *) 'Iteracoes: ', iter
      write(20, *) 'Convergiu: ', converged
   else
      write(20, *) 'Analise: LINEAR'
   end if

   write(20, *) '----------------------------------------'
   write(20, *) 'DESLOCAMENTOS (por no):'
   do i=1,nn
      write(20,'(A,I4,A,1PE12.4,A,1PE12.4,A,1PE12.4)') &
         'No ', i, ': ux=', disp(gdl(i,1)), '  uy=', disp(gdl(i,2)), &
         '  rz=', disp(gdl(i,3))
   end do

   write(20, *) '----------------------------------------'
   write(20, *) 'REACOES (apenas gdl prescritos):'
   do i=1,nbc
      write(20,'(A,I4,A,I2,A,1PE12.4)') 'No ', bcnode(i), &
         ' dof ', bcdof(i), '  R=', reaction(gdl(bcnode(i),bcdof(i)))
   end do

   write(20, *) '----------------------------------------'
   write(20, *) 'ESFORCOS LOCAIS ELASTICOS POR ELEMENTO:'
   write(20, *) 'Convencao: [N1,V1,M1,N2,V2,M2]'
   do e=1,ne
      n1=e1(e); n2=e2(e); matid=emat(e)
      emod=ea(matid); area=aa(matid); iner=ei(matid)
      length = sqrt((x(n2)-x(n1))**2 + (y(n2)-y(n1))**2)
      c = (x(n2)-x(n1))/length
      s = (y(n2)-y(n1))/length
      call frame_elastic_local(emod, area, iner, length, kl)
      call transform6(c, s, trans)
      g = [gdl(n1,1), gdl(n1,2), gdl(n1,3), gdl(n2,1), gdl(n2,2), gdl(n2,3)]
      ue = matmul(trans, disp(g))
      qel = matmul(kl, ue)
      write(20,'(A,I4,A,1PE12.4)') 'Elem ', e, ': Nax=', nax(e)
      write(20,'(6(1PE12.4,1X))') qel
   end do
   close(20)

contains

   integer function gdl(node_id, dof)
      integer, intent(in) :: node_id, dof
      gdl = 3*(node_id-1) + dof
   end function gdl

   subroutine build_stiffness(use_geo)
      logical, intent(in) :: use_geo
      integer :: row, col, elem, ni, nj, material, gdls(6)
      real(dp) :: le, ce, se, ke(6,6), kg(6,6), kt(6,6), keg(6,6), t(6,6)

      kmat = 0.0_dp
      do elem=1,ne
         ni=e1(elem); nj=e2(elem); material=emat(elem)
         le = sqrt((x(nj)-x(ni))**2 + (y(nj)-y(ni))**2)
         ce = (x(nj)-x(ni))/le
         se = (y(nj)-y(ni))/le
         call frame_elastic_local(ea(material), aa(material), ei(material), le, ke)
         if (use_geo) then
            call frame_geometric_local(nax(elem), le, kg)
         else
            kg = 0.0_dp
         end if
         kt = ke + kg
         call transform6(ce, se, t)
         keg = matmul(transpose(t), matmul(kt, t))
         gdls = [gdl(ni,1), gdl(ni,2), gdl(ni,3), gdl(nj,1), gdl(nj,2), gdl(nj,3)]
         do row=1,6
            do col=1,6
               kmat(gdls(row),gdls(col)) = kmat(gdls(row),gdls(col)) + keg(row,col)
            end do
         end do
      end do
   end subroutine build_stiffness

   subroutine frame_elastic_local(ee, area_e, iner_e, le, k)
      real(dp), intent(in) :: ee, area_e, iner_e, le
      real(dp), intent(out) :: k(6,6)
      real(dp) :: ea_l, ei_l, ei_l2, ei_l3

      k=0.0_dp
      ea_l = ee*area_e/le
      ei_l = ee*iner_e/le
      ei_l2 = ee*iner_e/le**2
      ei_l3 = ee*iner_e/le**3
      k(1,1)= ea_l; k(1,4)=-ea_l; k(4,1)=-ea_l; k(4,4)= ea_l
      k(2,2)= 12.0_dp*ei_l3; k(2,3)=  6.0_dp*ei_l2
      k(2,5)=-12.0_dp*ei_l3; k(2,6)=  6.0_dp*ei_l2
      k(3,2)=  6.0_dp*ei_l2; k(3,3)=  4.0_dp*ei_l
      k(3,5)= -6.0_dp*ei_l2; k(3,6)=  2.0_dp*ei_l
      k(5,2)=-12.0_dp*ei_l3; k(5,3)= -6.0_dp*ei_l2
      k(5,5)= 12.0_dp*ei_l3; k(5,6)= -6.0_dp*ei_l2
      k(6,2)=  6.0_dp*ei_l2; k(6,3)=  2.0_dp*ei_l
      k(6,5)= -6.0_dp*ei_l2; k(6,6)=  4.0_dp*ei_l
   end subroutine frame_elastic_local

   subroutine frame_geometric_local(axial, le, k)
      real(dp), intent(in) :: axial, le
      real(dp), intent(out) :: k(6,6)
      real(dp) :: fac

      k=0.0_dp
      fac = axial/(30.0_dp*le)
      k(2,2)= 36.0_dp*fac;    k(2,3)=  3.0_dp*le*fac
      k(2,5)=-36.0_dp*fac;    k(2,6)=  3.0_dp*le*fac
      k(3,2)=  3.0_dp*le*fac; k(3,3)=  4.0_dp*le**2*fac
      k(3,5)= -3.0_dp*le*fac; k(3,6)= -1.0_dp*le**2*fac
      k(5,2)=-36.0_dp*fac;    k(5,3)= -3.0_dp*le*fac
      k(5,5)= 36.0_dp*fac;    k(5,6)= -3.0_dp*le*fac
      k(6,2)=  3.0_dp*le*fac; k(6,3)= -1.0_dp*le**2*fac
      k(6,5)= -3.0_dp*le*fac; k(6,6)=  4.0_dp*le**2*fac
   end subroutine frame_geometric_local

   subroutine transform6(ce, se, t)
      real(dp), intent(in) :: ce, se
      real(dp), intent(out) :: t(6,6)
      t=0.0_dp
      t(1,1)= ce; t(1,2)= se; t(2,1)=-se; t(2,2)= ce; t(3,3)=1.0_dp
      t(4,4)= ce; t(4,5)= se; t(5,4)=-se; t(5,5)= ce; t(6,6)=1.0_dp
   end subroutine transform6

   subroutine solve_with_bc()
      integer :: ii, jj, gd, row, col
      real(dp) :: penalty

      disp = 0.0_dp
      do i=1,nbc
         disp(gdl(bcnode(i),bcdof(i))) = bcval(i)
      end do

      if (method == 1) then
         kwork = kmat
         fwork = force
         penalty = 0.0_dp
         do i=1,ndof
            penalty = max(penalty, abs(kwork(i,i)))
         end do
         if (penalty == 0.0_dp) penalty = 1.0_dp
         penalty = 1.0e12_dp*penalty
         do i=1,nbc
            gd = gdl(bcnode(i), bcdof(i))
            kwork(gd,gd) = kwork(gd,gd) + penalty
            fwork(gd) = fwork(gd) + penalty*bcval(i)
         end do
         call gauss(ndof, kwork, fwork, disp)
      else if (method == 2) then
         nfree = 0
         do i=1,ndof
            if (prescribed(i) == 0) then
               nfree = nfree + 1
               free_gdls(nfree) = i
            end if
         end do
         kwork = 0.0_dp
         fwork = 0.0_dp
         do ii=1,nfree
            row = free_gdls(ii)
            fwork(ii) = force(row)
            do col=1,ndof
               if (prescribed(col) == 1) fwork(ii) = fwork(ii) - kmat(row,col)*disp(col)
            end do
            do jj=1,nfree
               col = free_gdls(jj)
               kwork(ii,jj) = kmat(row,col)
            end do
         end do
         call gauss(nfree, kwork, fwork, ured)
         do ii=1,nfree
            disp(free_gdls(ii)) = ured(ii)
         end do
      else
         write(*,*) 'Erro: METHOD deve ser 1 ou 2.'
         stop
      end if
   end subroutine solve_with_bc

   subroutine compute_axials(u, axial)
      real(dp), intent(in) :: u(maxdof)
      real(dp), intent(out) :: axial(maxe)
      integer :: elem, ni, nj, material, gdls(6)
      real(dp) :: le, ce, se, t(6,6), ulocal(6)

      do elem=1,ne
         ni=e1(elem); nj=e2(elem); material=emat(elem)
         le = sqrt((x(nj)-x(ni))**2 + (y(nj)-y(ni))**2)
         ce = (x(nj)-x(ni))/le
         se = (y(nj)-y(ni))/le
         call transform6(ce, se, t)
         gdls = [gdl(ni,1), gdl(ni,2), gdl(ni,3), gdl(nj,1), gdl(nj,2), gdl(nj,3)]
         ulocal = matmul(t, u(gdls))
         axial(elem) = ea(material)*aa(material)/le*(ulocal(4)-ulocal(1))
      end do
   end subroutine compute_axials

   subroutine gauss(n, a, b, sol)
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(maxdof,maxdof), b(maxdof)
      real(dp), intent(out) :: sol(maxdof)
      integer :: k, row, col, pivot_row
      real(dp) :: pivot, factor, amax, tmp, sumv, tmp_row(maxdof)

      do k=1,n-1
         pivot_row = k
         amax = abs(a(k,k))
         do row=k+1,n
            if (abs(a(row,k)) > amax) then
               amax = abs(a(row,k))
               pivot_row = row
            end if
         end do
         if (amax < 1.0e-30_dp) then
            write(*,*) 'Pivo ~ 0 em K=', k, '; sistema singular.'
            stop
         end if
         if (pivot_row /= k) then
            tmp_row(1:n) = a(k,1:n)
            a(k,1:n) = a(pivot_row,1:n)
            a(pivot_row,1:n) = tmp_row(1:n)
            tmp = b(k); b(k) = b(pivot_row); b(pivot_row) = tmp
         end if
         pivot = a(k,k)
         do row=k+1,n
            factor = a(row,k)/pivot
            a(row,k) = 0.0_dp
            do col=k+1,n
               a(row,col) = a(row,col) - factor*a(k,col)
            end do
            b(row) = b(row) - factor*b(k)
         end do
      end do

      if (abs(a(n,n)) < 1.0e-30_dp) then
         write(*,*) 'Pivo final ~ 0; sistema singular.'
         stop
      end if
      sol(n) = b(n)/a(n,n)
      do row=n-1,1,-1
         sumv = 0.0_dp
         do col=row+1,n
            sumv = sumv + a(row,col)*sol(col)
         end do
         sol(row) = (b(row)-sumv)/a(row,row)
      end do
   end subroutine gauss

end program portico2d_pdelta
