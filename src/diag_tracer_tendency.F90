module mod_diag_tracer_tendency
   ! Total tracer tendency for the ltra_diag budget diagnostics [mmol/m2/s]:
   !    tra_tendency = (c*h at end of step - c*h at start of step) / dt
   ! store_tracer_old is called at the start of the step (before REcoM),
   ! tracer_tendency after oce_timestep_ale. Both use hnode, which holds the
   ! thickness that matches the tracer values at these points: update_thickness_ale
   ! sets hnode = hnode_new at the end of the step, while hnode_new does not hold
   ! the current thickness then, as with ldiag_DVD it gets the old thickness back.
   ! In single precision (use_tracer_carry) the layer content also includes the
   ! rounding carry: c*h + carry.
   USE o_PARAM, only: use_tracer_carry
   USE MOD_MESH
   USE MOD_TRACER
   USE MOD_PARTIT
   implicit none
   private
   public :: store_tracer_old, tracer_tendency, tracer_budget_ml

contains

subroutine store_tracer_old(tracers, partit, mesh)
   type(t_tracer), intent(inout) :: tracers
   type(t_partit), intent(in)    :: partit
   type(t_mesh)  , intent(in)    :: mesh

   integer :: n, nz, nu1, nl1, tr_num

   do tr_num=1, tracers%num_tracers
      if (.not. tracers%data(tr_num)%ltra_diag) cycle
      do n = 1, partit%myDim_nod2D
         nu1 = mesh%ulevels_nod2D(n)
         nl1 = mesh%nlevels_nod2D(n)
         do nz = nu1, nl1-1
            tracers%work%tracer_old(nz,n,tr_num) = tracers%data(tr_num)%values(nz,n)*mesh%hnode(nz,n)
         end do
         if (use_tracer_carry) tracers%work%tracer_old(nu1:nl1-1,n,tr_num) = &
              tracers%work%tracer_old(nu1:nl1-1,n,tr_num) + tracers%work%carry(nu1:nl1-1,n,tr_num)
      end do
   end do
end subroutine store_tracer_old

subroutine tracer_tendency(tracers, partit, mesh)
   use g_config, only : dt
   type(t_tracer), intent(inout) :: tracers
   type(t_partit), intent(in)    :: partit
   type(t_mesh)  , intent(in)    :: mesh

   integer :: n, nz, nu1, nl1, tr_num

   do tr_num=1, tracers%num_tracers
      if (.not. tracers%data(tr_num)%ltra_diag) cycle
      do n = 1, partit%myDim_nod2D
         nu1 = mesh%ulevels_nod2D(n)
         nl1 = mesh%nlevels_nod2D(n)
         do nz = nu1, nl1-1
            tracers%work%tra_tendency(nz,n,tr_num) = &
                 (tracers%data(tr_num)%values(nz,n)*mesh%hnode(nz,n) &
                - tracers%work%tracer_old(nz,n,tr_num)) / dt
         end do
         if (use_tracer_carry) tracers%work%tra_tendency(nu1:nl1-1,n,tr_num) = &
              tracers%work%tra_tendency(nu1:nl1-1,n,tr_num) + tracers%work%carry(nu1:nl1-1,n,tr_num)/dt
      end do
   end do
end subroutine tracer_tendency

subroutine tracer_budget_ml(tracers, partit, mesh)
   use o_ARRAYS, only : MLD3          
   type(t_tracer), intent(inout) :: tracers
   type(t_partit), intent(in)    :: partit
   type(t_mesh)  , intent(in)    :: mesh

   integer       :: n, nz, nu1, nl1, tr_num
   real(kind=WP) :: h, ztop, f, hml, tint
   real(kind=WP) :: s_tend, s_hadv, s_vadv, s_hdiff, s_vdiff, s_sms

   do tr_num=1, tracers%num_tracers
      if (.not. tracers%data(tr_num)%ltra_diag) cycle
      do n = 1, partit%myDim_nod2D
         nu1 = mesh%ulevels_nod2D(n)
         nl1 = mesh%nlevels_nod2D(n)
         ! MLD3 is negative downward
         ! Smallest ML = top layer
         h    = max(abs(MLD3(n)), mesh%hnode(nu1,n))
         ztop = 0.0_WP
         hml  = 0.0_WP; tint = 0.0_WP
         s_tend=0.0_WP; s_hadv=0.0_WP; s_vadv=0.0_WP; s_hdiff=0.0_WP; s_vdiff=0.0_WP; s_sms=0.0_WP
         do nz = nu1, nl1-1
            f = max(0.0_WP, min(1.0_WP, (h-ztop)/mesh%hnode(nz,n)))
            if (f <= 0.0_WP) exit
            hml   = hml  + f*mesh%hnode(nz,n)
            tint  = tint + f*mesh%hnode(nz,n)*tracers%data(tr_num)%values(nz,n)
            s_tend  = s_tend  + f*tracers%work%tra_tendency(nz,n,tr_num)
            s_hadv  = s_hadv  + f*tracers%work%tra_advhoriz(nz,n,tr_num)
            s_vadv  = s_vadv  + f*tracers%work%tra_advvert(nz,n,tr_num)
            s_hdiff = s_hdiff + f*( tracers%work%tra_diff_part_hor_redi(nz,n,tr_num) &
                                  + tracers%work%tra_diff_part_ver_redi_expl(nz,n,tr_num) )
            if (tracers%data(tr_num)%i_vert_diff) then
               s_vdiff = s_vdiff + f*tracers%work%tra_diff_part_ver_impl(nz,n,tr_num)
            else
               s_vdiff = s_vdiff + f*tracers%work%tra_diff_part_ver_expl(nz,n,tr_num)
            end if
            ztop = ztop + mesh%hnode(nz,n)
         end do
         if (hml > 0.0_WP) then
            tracers%work%ml_h(n) = hml
            tracers%work%ml_tmean(n,tr_num) = tint/hml
         
            tracers%work%ml_tend(n,tr_num)     = s_tend/hml
            tracers%work%ml_hor_adv(n,tr_num)  = s_hadv/hml
            tracers%work%ml_ver_adv(n,tr_num)  = s_vadv/hml
            tracers%work%ml_redi_diff(n,tr_num) = s_hdiff/hml
            tracers%work%ml_ver_diff(n,tr_num) = s_vdiff/hml
            tracers%work%ml_resid(n,tr_num)    = (s_tend - (s_hadv + s_vadv + s_hdiff + s_vdiff + s_sms))/hml

         else
                tracers%work%ml_h(n) = 0.0_WP
                tracers%work%ml_tmean(n,tr_num) = 0.0_WP
                tracers%work%ml_tend(n,tr_num)  = 0.0_WP
                tracers%work%ml_hor_adv(n,tr_num) = 0.0_WP
                tracers%work%ml_ver_adv(n,tr_num) = 0.0_WP
                tracers%work%ml_redi_diff(n,tr_num) = 0.0_WP
                tracers%work%ml_ver_diff(n,tr_num) = 0.0_WP
                tracers%work%ml_resid(n,tr_num)    = 0.0_WP
         end if
      end do
   end do
end subroutine tracer_budget_ml

end module mod_diag_tracer_tendency
