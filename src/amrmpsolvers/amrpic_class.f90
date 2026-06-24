!> AMRPIC solver class
!> Provides piecewise interface calculation (PIC) type for AMRVOF class
!> Picks between IRL-free or IRL-based implementation
#ifdef USE_IRL
module amrpic_class
! IRL: use IRL's separator union multifab to define the interface
   use amrex_amr_module,          only: WP => amrex_real 
   use amrex_sepunionmfab_module, only: picmfab_build => amrex_sepunionmfab_build, &
                                        picmfab_destroy => amrex_sepunionmfab_destroy, &
                                        picmfab_rebuild => amrex_sepunionmfab_rebuild, &
                                        picmfab => amrex_sepunionmfab
    use irl_fortran_interface,    only: pic_container => SeparatorUnion_type_raw

   implicit none
   
   integer :: pic_ncomp = 1

   interface reflect_pic
      module procedure reflect_pic_raw_doubles
      module procedure reflect_pic_sepunion
   end interface

   interface get_plane
      module procedure get_plane_raw_doubles
      module procedure get_plane_sepunion
   end interface

   interface set_to_one_plane
      module procedure set_to_one_plane_raw_doubles
      module procedure set_to_one_plane_sepunion
   end interface

   interface set_to_full
      module procedure set_to_full_raw_doubles
      module procedure set_to_full_sepunion
   end interface

   interface set_to_empty
      module procedure set_to_empty_raw_doubles
      module procedure set_to_empty_sepunion
   end interface

   interface set_to_full_or_empty
      module procedure set_to_full_or_empty_raw_doubles
      module procedure set_to_full_or_empty_sepunion
   end interface

   interface is_full
      module procedure is_full_raw_doubles
      module procedure is_full_sepunion
   end interface

   interface is_empty
      module procedure is_empty_raw_doubles
      module procedure is_empty_sepunion
   end interface

   interface is_one_plane
      module procedure is_one_plane_raw_doubles
      module procedure is_one_plane_sepunion
   end interface

   interface is_paraboloid
      module procedure is_paraboloid_raw_doubles
      module procedure is_paraboloid_sepunion
   end interface

   interface cut_rectcub_pic
      module procedure cut_rectcub_pic_raw_doubles
      module procedure cut_rectcub_pic_sepunion
   end interface

   interface cut_tet_pic
      module procedure cut_tet_pic_raw_doubles
      module procedure cut_tet_pic_sepunion
   end interface

   contains

   function get_normal(plane) result(normal)
        implicit none
        real(WP), dimension(4), intent(in)  :: plane
        real(WP), dimension(3) :: normal
        normal = plane(1:3)
   end function get_normal

   subroutine reflect_pic_raw_doubles(this, ref_pic, dir, loc)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      real(WP), dimension(:), intent(in) :: ref_pic
      integer, intent(in) :: dir
      real(WP), intent(in) :: loc
      this=ref_pic
      this(dir)=-this(dir)
      this(4)=this(4)-2.0_WP*ref_pic(dir)*loc
   end subroutine reflect_pic_raw_doubles

   subroutine reflect_pic_sepunion(this, ref_pic, dir, loc)
      use irl_fortran_interface, only: reflect
      implicit none
      type(pic_container), intent(inout) :: this(..)
      type(pic_container), intent(in) :: ref_pic(..)
      integer, intent(in) :: dir
      real(WP), intent(in) :: loc
      integer :: i 
      select rank(this)
         rank(0)
            select rank(ref_pic)
               rank(0)
                  call reflect(this, ref_pic, dir, loc)
               rank(1)
                  call reflect(this, ref_pic(lbound(ref_pic,1)), dir, loc)
            end select
         rank(1)
            select rank(ref_pic)
               rank(0)
                  do i = lbound(this,1), ubound(this,1)
                     call reflect(this(i), ref_pic, dir, loc)
                  end do
               rank(1)
                  do i = max(lbound(this,1),lbound(ref_pic,1)), min(ubound(this,1),ubound(ref_pic,1))
                     call reflect(this(i), ref_pic(i), dir, loc)
                  end do
            end select
      end select
   end subroutine reflect_pic_sepunion

   function get_plane_raw_doubles(this) result(a_plane_listed)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      real(WP), dimension(4) :: a_plane_listed
      a_plane_listed = this
   end function get_plane_raw_doubles

   function get_plane_sepunion(this) result(a_plane_listed)
      use irl_fortran_interface
      implicit none
      type(pic_container), intent(inout) :: this(..)
      real(WP), dimension(4) :: a_plane_listed
      select rank(this)
         rank(0)
            a_plane_listed = getPlane(this,0)
         rank(1)
            a_plane_listed = getPlane(this(lbound(this,1)),0)
      end select
   end function get_plane_sepunion

   subroutine set_to_one_plane_raw_doubles(this, plane_normal, plane_distance)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      real(WP), dimension(1:3), intent(in) :: plane_normal
      real(WP), intent(in) :: plane_distance
      this(1:3) = plane_normal
      this(4) = plane_distance
   end subroutine set_to_one_plane_raw_doubles

   subroutine set_to_one_plane_sepunion(this, plane_normal, plane_distance)
      use irl_fortran_interface, only: setToOnePlane
      implicit none
      type(pic_container), intent(inout) :: this(..)
      real(WP), dimension(1:3), intent(in) :: plane_normal
      real(WP), intent(in) :: plane_distance
      integer :: i
      select rank(this)
         rank(0)
            call setToOnePlane(this,plane_normal,plane_distance)
         rank(1)
            do i = lbound(this,1), ubound(this,1)
               call setToOnePlane(this(i),plane_normal,plane_distance)
            end do 
      end select
   end subroutine set_to_one_plane_sepunion

   subroutine set_to_full_raw_doubles(this)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      this(1:3) = [0.0_WP,0.0_WP,0.0_WP]
      this(4) = 1.0e15_WP
   end subroutine set_to_full_raw_doubles 

   subroutine set_to_full_sepunion(this)
      use irl_fortran_interface, only: setToFull
      implicit none
      type(pic_container), intent(inout) :: this(..)
      integer :: i,j,k,l 
      select rank(this)
         rank(0)
            call setToFull(this)
         rank(1)
            do i = lbound(this,1), ubound(this,1)
               call setToFull(this(i))
            end do
         rank(4)
            do i = lbound(this,1), ubound(this,1)
               do j = lbound(this,2), ubound(this,2)
                  do k = lbound(this,3), ubound(this,3)
                     do l = lbound(this,4), ubound(this,4)
                        call setToFull(this(i,j,k,l))
                     end do
                  end do
               end do
            end do
      end select
   end subroutine set_to_full_sepunion 

   subroutine set_to_empty_raw_doubles(this)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      this(1:3) = [0.0_WP,0.0_WP,0.0_WP]
      this(4) = -1.0e15_WP
   end subroutine set_to_empty_raw_doubles 

   subroutine set_to_empty_sepunion(this)
      use irl_fortran_interface, only: setToEmpty
      implicit none
      type(pic_container), intent(inout) :: this(..)
      integer :: i,j,k,l
      select rank(this)
         rank(0)
            call setToEmpty(this)
         rank(1)
            do i = lbound(this,1), ubound(this,1)
               call setToEmpty(this(i))
            end do
         rank(4)
            do i = lbound(this,1), ubound(this,1)
               do j = lbound(this,2), ubound(this,2)
                  do k = lbound(this,3), ubound(this,3)
                     do l = lbound(this,4), ubound(this,4)
                        call setToEmpty(this(i,j,k,l))
                     end do
                  end do 
               end do
            end do 
      end select
   end subroutine set_to_empty_sepunion 

   subroutine set_to_full_or_empty_raw_doubles(this, VF)
      implicit none
      real(WP), dimension(:), intent(inout) :: this
      real(WP), intent(in) :: VF
      this(1:3) = [0.0_WP,0.0_WP,0.0_WP]
      this(4) = sign(1.0e15_WP,VF-0.5_WP)
   end subroutine set_to_full_or_empty_raw_doubles

   subroutine set_to_full_or_empty_sepunion(this, VF)
      implicit none
      type(pic_container), intent(inout) :: this(..)
      real(WP), intent(in) :: VF
      integer :: i,j,k,l 
      if (VF.gt.0.5_WP) then 
         call set_to_full(this)
      else 
         call set_to_empty(this) 
      end if 
   end subroutine set_to_full_or_empty_sepunion

   function is_full_raw_doubles(this) result(res)
      implicit none
      real(WP), dimension(..), intent(in) :: this
      logical :: res
      integer :: i,j,k
      res = .true.
      select rank(this)
         rank(0)
            if (this.le.+1.0e9_WP) res = .false.
         rank(1)
            if (this(4).le.+1.0e9_WP) res = .false.
         rank(4)
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     if (this(i,j,k,4).le.+1.0e9_WP) res = .false.
                  end do 
               end do 
            end do
      end select
   end function is_full_raw_doubles

   function is_full_sepunion(this) result(res)
      use irl_fortran_interface, only: isFull
      implicit none
      type(pic_container), intent(inout) :: this(..)
      logical :: res
      integer :: i,j,k,l
      res = .true.
      select rank(this)
         rank(0)
            res = isFull(this)
         rank(1)
            res = isFull(this(lbound(this,1)))
         rank(4)
            res = .true.
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     do l=lbound(this,4),ubound(this,4)
                        res = res.and.isFull(this(i,j,k,l))
                     end do 
                  end do 
               end do 
            end do
      end select
   end function is_full_sepunion

   function is_empty_raw_doubles(this) result(res)
      implicit none
      real(WP), dimension(..), intent(in) :: this
      logical :: res
      integer :: i,j,k
      res = .true.
      select rank(this)
         rank(0)
            if (this.ge.-1.0e9_WP) res = .false.
         rank(1)
            if (this(4).ge.-1.0e9_WP) res = .false.
         rank(4)
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     if (this(i,j,k,4).ge.-1.0e9_WP) res = .false.
                  end do 
               end do 
            end do
      end select
   end function is_empty_raw_doubles

   function is_empty_sepunion(this) result(res)
      use irl_fortran_interface, only: isEmpty
      implicit none
      type(pic_container), intent(inout) :: this(..)
      logical :: res
      integer :: i,j,k,l
      res = .true.
      select rank(this)
         rank(0)
            res = isEmpty(this)
         rank(1)
            res = isEmpty(this(lbound(this,1)))
         rank(4)
            res = .true.
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     do l=lbound(this,4),ubound(this,4)
                        res = res.and.isEmpty(this(i,j,k,l))
                     end do 
                  end do 
               end do 
            end do
      end select
   end function is_empty_sepunion

   function is_one_plane_raw_doubles(this) result(res)
      implicit none
      real(WP), dimension(:), intent(in) :: this
      logical :: res
      res = .true.
   end function is_one_plane_raw_doubles

   function is_one_plane_sepunion(this) result(res)
      use irl_fortran_interface, only: isOnePlane
      implicit none
      type(pic_container), intent(inout) :: this(..)
      logical :: res
      integer :: i,j,k,l
      res = .true.
      select rank(this)
         rank(0)
            res = isOnePlane(this)
         rank(1)
            res = isOnePlane(this(lbound(this,1)))
         rank(4)
            res = .true.
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     do l=lbound(this,4),ubound(this,4)
                        res = res.and.isOnePlane(this(i,j,k,l))
                     end do 
                  end do 
               end do 
            end do
      end select
   end function is_one_plane_sepunion

   function is_paraboloid_raw_doubles(this) result(res)
      implicit none
      real(WP), dimension(:), intent(in) :: this
      logical :: res
      res = .false.
   end function is_paraboloid_raw_doubles

   function is_paraboloid_sepunion(this) result(res)
      use irl_fortran_interface, only: isParaboloid
      implicit none
      type(pic_container), intent(inout) :: this(..)
      logical :: res
      integer :: i,j,k,l
      res = .true.
      select rank(this)
         rank(0)
            res = isParaboloid(this)
         rank(1)
            res = isParaboloid(this(lbound(this,1)))
         rank(4)
            res = .true.
            do i=lbound(this,1),ubound(this,1)
               do j=lbound(this,2),ubound(this,2)
                  do k=lbound(this,3),ubound(this,3)
                     do l=lbound(this,4),ubound(this,4)
                        res = res.and.isParaboloid(this(i,j,k,l))
                     end do 
                  end do 
               end do 
            end do
      end select
   end function is_paraboloid_sepunion

   subroutine cut_rectcub_pic_raw_doubles(ptlo, pthi, pic, vol_liq, vol_gas, bary_liq, bary_gas)
      use amrvof_geometry, only: cut_hex_vol
      implicit none
      real(WP), dimension(3), intent(in) :: ptlo, pthi
      real(WP), dimension(:), intent(in) :: pic
      real(WP), intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3), intent(out) :: bary_liq, bary_gas
      real(WP), dimension(3,8) :: hex
      ! Build hex cell
      hex(:,1)=[ptlo(1),ptlo(2),ptlo(3)]; hex(:,2)=[pthi(1),ptlo(2),ptlo(3)] 
      hex(:,3)=[pthi(1),pthi(2),ptlo(3)]; hex(:,4)=[ptlo(1),pthi(2),ptlo(3)]
      hex(:,5)=[ptlo(1),ptlo(2),pthi(3)]; hex(:,6)=[pthi(1),ptlo(2),pthi(3)]
      hex(:,7)=[pthi(1),pthi(2),pthi(3)]; hex(:,8)=[ptlo(1),pthi(2),pthi(3)]
      ! Cut hex by plane
      call cut_hex_vol(hex,pic,vol_liq,vol_gas,bary_liq,bary_gas)
   end subroutine cut_rectcub_pic_raw_doubles

   subroutine cut_rectcub_pic_sepunion(ptlo, pthi, pic, vol_liq, vol_gas, bary_liq, bary_gas)
      use irl_fortran_interface, only: RectCub_type,SepVM_type,new,construct_2pt,getNormMoments,getVolume,getCentroid
      implicit none
      real(WP), dimension(3), intent(in) :: ptlo, pthi
      type(pic_container), intent(in) :: pic(..)
      real(WP), intent(out) :: vol_liq
      real(WP), intent(out), optional :: vol_gas
      real(WP), dimension(3), intent(out), optional :: bary_liq
      real(WP), dimension(3), intent(out), optional :: bary_gas
      type(RectCub_type) :: cell
      type(SepVM_type) :: vol_and_bary
      call new(cell)
      call new(vol_and_bary)
      ! Construct rectangular cuboid from points
      call construct_2pt(cell,ptlo,pthi)
      ! Compute normalized separated volume-moments using IRL
      select rank(pic)
         rank(0)     
            call getNormMoments(cell,pic,vol_and_bary)
         rank(1) ! Get VM associated with first element in array of PICs              
            call getNormMoments(cell,pic(lbound(pic,1)),vol_and_bary)
      end select
      ! Extract moments from IRL object
      vol_liq  = getVolume  (vol_and_bary,0)
      if (present(vol_gas))  vol_gas  = getVolume  (vol_and_bary,1)
      if (present(bary_liq)) bary_liq = getCentroid(vol_and_bary,0)
      if (present(bary_gas)) bary_gas = getCentroid(vol_and_bary,1)
   end subroutine cut_rectcub_pic_sepunion

   subroutine cut_tet_pic_raw_doubles(tet, pic, VF0, vol_tot, bary_tot, flux)
      use amrvof_geometry, only: cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert,cut_nntet,tet_vol
      implicit none
      real(WP), dimension(3,4), intent(in) :: tet
      real(WP), dimension(:), intent(in) :: pic
      real(WP), intent(in) :: VF0,vol_tot
      real(WP), dimension(3), intent(in) :: bary_tot
      real(WP), dimension(8), intent(inout) :: flux
      real(WP), dimension(3) :: a,b,c,bary,normal
      real(WP) :: mu,my_vol,dist
      real(WP), dimension(4) :: dd
      real(WP), dimension(3,8) :: vert
      integer :: icase,n1,v1,v2
      ! Get PLIC from this cell
      normal=pic(1:3)
      dist  =pic( 4 )
      
      ! Compute signed distance to plane for each vertex
      dd(1)=normal(1)*tet(1,1)+normal(2)*tet(2,1)+normal(3)*tet(3,1)-dist
      dd(2)=normal(1)*tet(1,2)+normal(2)*tet(2,2)+normal(3)*tet(3,2)-dist
      dd(3)=normal(1)*tet(1,3)+normal(2)*tet(2,3)+normal(3)*tet(3,3)-dist
      dd(4)=normal(1)*tet(1,4)+normal(2)*tet(2,4)+normal(3)*tet(3,4)-dist
      
      ! Find cut case
      icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
      &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
      &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
      &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))
      
      ! Copy vertices
      vert(:,1:4)=tet(:,1:4)
      
      ! Create interpolated vertices on cut plane
      do n1=1,cut_nvert(icase)
         v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
         mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
         vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
      end do
      
      ! Cut the minority phase (safer as we subtract small from large)
      if (VF0.gt.0.5_WP) then
         ! Liquid is dominant → compute gas directly
         do n1=1,cut_nntet(icase)-1
            a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
            &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
            flux( 2 )=flux( 2 )+my_vol
            flux(6:8)=flux(6:8)+my_vol*bary
         end do
         ! Liquid=total-gas
         flux( 1 )=         vol_tot-flux( 2 )
         flux(3:5)=bary_tot*vol_tot-flux(6:8)
      else
         ! Gas is dominant → compute liquid directly
         do n1=cut_ntets(icase),cut_nntet(icase),-1
            a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
            &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
            flux( 1 )=flux( 1 )+my_vol
            flux(3:5)=flux(3:5)+my_vol*bary
         end do
         ! Gas=total-liquid
         flux( 2 )=         vol_tot-flux( 1 )
         flux(6:8)=bary_tot*vol_tot-flux(3:5)
      end if
   end subroutine cut_tet_pic_raw_doubles

   subroutine cut_tet_pic_sepunion(tet, pic, VF0, vol_tot, bary_tot, flux)
      use irl_fortran_interface, only: Tet_type,SepVM_type,new,construct,getNormMoments,getVolume,getCentroid,getPlane
      use amrvof_geometry, only: cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert,cut_nntet,tet_vol
      implicit none
      real(WP), dimension(3,4), intent(in) :: tet
      type(pic_container), intent(in) :: pic(..)
      real(WP), intent(in) :: VF0,vol_tot
      real(WP), dimension(3), intent(in) :: bary_tot
      real(WP), dimension(8), intent(inout) :: flux
      type(Tet_type) :: cell
      type(SepVM_type) :: vol_and_bary
      call new(cell)
      call new(vol_and_bary)
      ! Construct rectangular cuboid from points
      call construct(cell,tet)
      ! Compute normalized separated volume-moments using IRL
      select rank(pic)
         rank(0)                   
            call getNormMoments(cell,pic,vol_and_bary)
         rank(1) ! Get VM associated with first element in array of PICs              
            call getNormMoments(cell,pic(lbound(pic,1)),vol_and_bary)
      end select
      ! Extract moments from IRL object
      flux( 1 ) = abs(getVolume(vol_and_bary,0))
      flux( 2 ) = abs(getVolume(vol_and_bary,1))
      flux(3:5) = flux(1)*getCentroid(vol_and_bary,0)
      flux(6:8) = flux(2)*getCentroid(vol_and_bary,1)
   end subroutine cut_tet_pic_sepunion

end module amrpic_class
#endif
