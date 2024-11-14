! flow_field.f90
!
!****************************************************************************
!
!   This file is part of the Yardang Dynamics Model.
!
!   Calculate the flow field across the topography.
!
!   @Copyright Key Laboratory of Earth System Numerical Modeling and Application, 
!   University of Chinese Academy of Sciences, Beijing, China.
!   All rights reserved.
!
!   Created by Haoxuan Dang, 2024.
!
!   When using the code, please cite the associated papers.
!
!****************************************************************************
    

module flow_field
    
    implicit none
    
    contains
    
    subroutine set_zero(ia,ja,h,u,v,p,q,qu,uu,pv,vv,hu,hv)

        implicit none
    
        integer ia,ja
        real(kind=8) h(ia+2,ja+2,2),hu(ia+3,ja+2,2),hv(ia+2,ja+3,2)
        real(kind=8) u(ia+3,ja+2,2),p(ia+3,ja+2,2),v(ia+2,ja+3,2),q(ia+2,ja+3,2)
        real(kind=8) qu(ia+1,ja+1,2),uu(ia+1,ja+1,2),pv(ia+1,ja+1,2),vv(ia+1,ja+1,2)
    
        h = 0.d0; u = 0.d0; v = 0.d0; p = 0.d0; q = 0.d0; qu = 0.d0
        uu = 0.d0; pv = 0.d0; vv = 0.d0; hu = 0.d0; hv = 0.d0
    
    end subroutine set_zero
    
    
    subroutine cal_hu(ia,ja,m,n,tc,td,h_min,h,u,hu)       ! calculate the water depth at node u
    
        implicit none
    
        integer i,j,ia,ja,m,n,tc,td,mm,nn
        real(kind=8) h_min,rph(m+1,n,1),h(ia+2,ja+2,2),u(ia+3,ja+2,2),phi(m+1,n,1),hu(ia+3,ja+2,2)
    
        ! TVD scheme
        call cal_rph(ia,ja,m,n,tc,td,u,h,rph)
    
        mm = size(rph,1)
        nn = size(rph,2)
        call lim_fun_Q(mm,nn,rph,phi)
    
        do i = 2,ia+2
            do j = 2,ja+1     
                if (u(i,j,td) >= 0.d0) then
                    hu(i,j,tc) = h(i-1,j,tc) + 0.5d0 * phi(i,j,1) * (h(i,j,tc) - h(i-1,j,tc))
                else
                    hu(i,j,tc) = h(i,j,tc) + 0.5d0 * phi(i,j,1) * (h(i-1,j,tc) - h(i,j,tc))
                end if
            end do
        end do

        ! free flow boundary
        hu(2:ia+2,1,tc) = hu(2:ia+2,2,tc)
        hu(2:ia+2,ja+2,tc) = hu(2:ia+2,ja+1,tc)
        hu(1,1:ja+2,tc) = hu(2,1:ja+2,tc)
        hu(ia+3,1:ja+2,tc) = hu(ia+2,1:ja+2,tc)
    
        do i = 1,ia+3
            do j = 1,ja+2
                if (hu(i,j,tc) < h_min) then
                    hu(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_hu
    
    
    subroutine cal_huc(ia,ja,tc,h_min,h,huc)       
    
        implicit none
    
        integer i,j,ia,ja,tc
        real(kind=8) h_min,h(ia+2,ja+2,2),huc(ia+3,ja+2,2)
    
        ! CD
        do i = 2,ia+2
            do j = 2,ja+1     
                huc(i,j,tc) = (h(i-1,j,tc) + h(i,j,tc)) / 2d0
            end do
        end do

        ! free flow boundary
        huc(2:ia+2,1,tc) = huc(2:ia+2,2,tc)
        huc(2:ia+2,ja+2,tc) = huc(2:ia+2,ja+1,tc)
        huc(1,1:ja+2,tc) = huc(2,1:ja+2,tc)
        huc(ia+3,1:ja+2,tc) = huc(ia+2,1:ja+2,tc)
    
        do i = 1,ia+3
            do j = 1,ja+2
                if (huc(i,j,tc) < h_min) then
                    huc(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_huc
    
    
    subroutine cal_hv(ia,ja,m,n,tc,td,h_min,h,v,hv)     ! calculate the water depth at node v
    
        implicit none
    
        integer i,j,ia,ja,m,n,tc,td,mm,nn
        real(kind=8) h_min,rqh(m,n+1,1),h(ia+2,ja+2,2),v(ia+2,ja+3,2),phi(m,n+1,1),hv(ia+2,ja+3,2)
    
        ! TVD scheme
        call cal_rqh(ia,ja,m,n,tc,td,v,h,rqh)
    
        mm = size(rqh,1)
        nn = size(rqh,2)
        call lim_fun_Q(mm,nn,rqh,phi)
    
        do i=2,ia+1
            do j = 2,ja+2
                if (v(i,j,td) >= 0.d0) then
                    hv(i,j,tc) = h(i,j-1,tc) + 0.5d0 * phi(i,j,1) * (h(i,j,tc) - h(i,j-1,tc))
                else
                    hv(i,j,tc) = h(i,j,tc) + 0.5d0 * phi(i,j,1) * (h(i,j-1,tc) - h(i,j,tc))
                end if
            end do
        end do

        ! free flow boundary
        hv(1,2:ja+2,tc) = hv(2,2:ja+2,tc)
        hv(ia+2,2:ja+2,tc) = hv(ia+1,2:ja+2,tc)
        hv(1:ia+2,1,tc) = hv(1:ia+2,2,tc)
        hv(1:ia+2,ja+3,tc) = hv(1:ia+2,ja+2,tc)
    
        do i = 1,ia+2
            do j = 1,ja+3
                if (hv(i,j,tc) < h_min) then
                    hv(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_hv
   
    
    subroutine cal_hvc(ia,ja,tc,h_min,h,hvc)     
    
        implicit none
    
        integer i,j,ia,ja,tc
        real(kind=8) h_min,h(ia+2,ja+2,2),hvc(ia+2,ja+3,2)
    
        ! CD
        do i=2,ia+1
            do j = 2,ja+2
                hvc(i,j,tc) = (h(i,j-1,tc) + h(i,j,tc)) / 2d0
            end do
        end do

        ! free flow boundary
        hvc(1,2:ja+2,tc) = hvc(2,2:ja+2,tc)
        hvc(ia+2,2:ja+2,tc) = hvc(ia+1,2:ja+2,tc)
        hvc(1:ia+2,1,tc) = hvc(1:ia+2,2,tc)
        hvc(1:ia+2,ja+3,tc) = hvc(1:ia+2,ja+2,tc)
    
        do i = 1,ia+2
            do j = 1,ja+3
                if (hvc(i,j,tc) < h_min) then
                    hvc(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_hvc
    
    
    subroutine cal_p(ia,ja,p_ini,u,hu,p)       ! calculate the flux p

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) p_ini,u(ia+3,ja+2,2),hu(ia+3,ja+2,2),p(ia+3,ja+2,2)
    
        do i = 1,ia+3
            do j = 1,ja+2
                if (i == 1) then
                    p(i,j,1) = p_ini
                else
                    p(i,j,1) = hu(i,j,1) * u(i,j,1)
                end if
            end do
        end do
    
    end subroutine cal_p


    subroutine cal_q(ia,ja,q_ini,v,hv,q)       ! calculate the flux q

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) q_ini,v(ia+2,ja+3,2),hv(ia+2,ja+3,2),q(ia+2,ja+3,2)
    
        do i = 1,ia+2
            do j = 1,ja+3
                if (i == 1) then
                    q(i,j,1) = q_ini
                else
                    q(i,j,1) = hv(i,j,1) * v(i,j,1)
                end if
            end do
        end do
    
    end subroutine cal_q    
    
    
    subroutine cal_dt_and_C(ia,ja,dx,dy,g,dtmin,t,u,v,p,q,h,dt,C)

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) w,e,s,n,dtmin,vel,dt1,dt,t,g,C,dx,dy,nu1,nu2
        real(kind=8) u(ia+3,ja+2,2),p(ia+3,ja+2,2),v(ia+2,ja+3,2),q(ia+2,ja+3,2),h(ia+2,ja+2,2),deltat(ia,ja,1),nu(ia,ja,1)
    
        ! calculate the time step
        do i = 1,ia
            do j = 1,ja
                if ((-1d0) * u(i+1,j+1,1) > 0.d0) then
                    w = (-1d0) * u(i+1,j+1,1)
                else
                    w = 0.d0
                end if
                if (u(i+2,j+1,1) > 0.d0) then
                    e = u(i+2,j+1,1)
                else
                    e = 0.d0
                end if
                if ((-1d0) * v(i+1,j+1,1) > 0.d0) then
                    s = (-1d0) * v(i+1,j+1,1)
                else
                    s = 0.d0
                end if
                if (v(i+1,j+2,1) > 0.d0) then
                    n = v(i+1,j+2,1)
                else
                    n = 0.d0
                end if
                vel = w + e + s + n
                if (vel == 0.d0) then
                    deltat(i,j,1) = dtmin
                else
                    deltat(i,j,1) = dx / vel
                end if
            end do
        end do
        dt1 = minval(deltat)
        dt = min(dtmin,dt1)

        ! calculate the Courant number    
        do i = 1,ia
            do j = 1,ja
                if (h(i+1,j+1,1) == 0.d0) then
                    nu(i,j,1) = 0.d0
                else
                    nu1 = abs(q(i+1,j+2,1) + q(i+1,j+1,1)) / (2d0 * h(i+1,j+1,1)) + sqrt(g * h(i+1,j+1,1))
                    nu2 = abs(p(i+2,j+1,1) + p(i+1,j+1,1)) / (2d0 * h(i+1,j+1,1)) + sqrt(g * h(i+1,j+1,1))
                    nu(i,j,1) = max(nu1,nu2)
                end if
            end do
        end do
    
        C = dt / dx * maxval(nu)
    
        if (C > 1d0) then
            print *,"calculations begin to destabilize at t=",t
            stop
        end if
    
    end subroutine cal_dt_and_C
    
    
    subroutine cal_h(ia,ja,h_min,dt,dx,dy,h_in,h_out,h,p,q)

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dt,dx,dy,h_min
        real(kind=8) h_in(ja+2,1),h_out(ja+2,1),h(ia+2,ja+2,2),p(ia+3,ja+2,2),q(ia+2,ja+3,2)

        do i = 1,ia+2
            do j = 1,ja+2
                h(i,j,2) = h(i,j,1) - dt / dx * (p(i+1,j,1) - p(i,j,1)) - &
                    dt / dy * (q(i,j+1,1) - q(i,j,1))
            end do
        end do
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (h(i,j,2) < h_min) then
                    h(i,j,2) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_h
    
    
    subroutine cal_eta(ia,ja,h,B,eta)           ! calculate the free surface for the new time step

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) h(ia+2,ja+2,2),B(ia+2,ja+2,2),eta(ia+2,ja+2,1)

        do i = 1,ia+2
            do j = 1,ja+2
                eta(i,j,1) = h(i,j,2) + B(i,j,1)
            end do
        end do
 
    end subroutine cal_eta

    
    subroutine cal_uh_and_ph(ia,ja,m,n,tc,td,h_min,h,u,p,uh,ph)         ! calculate the flux p and velocity u at node h

        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) h(ia+2,ja+2,2),u(ia+3,ja+2,2),p(ia+3,ja+2,2),phd(ia+2,ja+2,2),uh(ia+2,ja+2,1),ph(ia+2,ja+2,1)
        real(kind=8) h_min,ru(m-1,n,1),phi(m-1,n,1)

        ! determine the flow direction
        do i = 1,ia+2
            do j = 1,ja+2
                phd(i,j,td) = (p(i+1,j,td) + p(i,j,td)) / 2d0
            end do
        end do
    
        call cal_rp(ia,ja,m,n,tc,td,phd,u,ru)
    
        mm = size(ru,1)
        nn = size(ru,2)
        call lim_fun_U(mm,nn,ru,phi)
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (phd(i,j,td) >= 0.d0) then
                    uh(i,j,tc) = u(i,j,tc) + 0.5d0 * phi(i,j,1) * (u(i+1,j,tc) - u(i,j,tc))
                else
                    uh(i,j,tc) = u(i+1,j,tc) + 0.5d0 * phi(i,j,1) * (u(i,j,tc) - u(i+1,j,tc))
                end if
            end do
        end do

        do i = 1,ia+2
            do j = 1,ja+2
                ph(i,j,tc) = (p(i+1,j,tc) + p(i,j,tc)) / 2d0
            end do
        end do

        do i = 1,ia+2
            do j = 1,ja+2
                if (h(i,j,tc) < h_min) then
                    ph(i,j,tc) = 0.d0
                    uh(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_uh_and_ph        

    
    subroutine cal_vh_and_qh(ia,ja,m,n,tc,td,h_min,h,v,q,vh,qh)            ! calculate the flux q and velocity v at node h

        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) h(ia+2,ja+2,2),v(ia+2,ja+3,2),q(ia+2,ja+3,2),qhd(ia+2,ja+2,2),vh(ia+2,ja+2,1),qh(ia+2,ja+2,1)
        real(kind=8) h_min,rv(m,n-1,1),phi(m,n-1,1)

        ! determine the flow direction
        do i = 1,ia+2
            do j = 1,ja+2
                qhd(i,j,td) = (q(i,j,td) + q(i,j+1,td)) / 2d0
            end do
        end do
    
        call cal_rq(ia,ja,m,n,tc,td,qhd,v,rv)
    
        mm = size(rv,1)
        nn = size(rv,2)
        call lim_fun_U(mm,nn,rv,phi)
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (qhd(i,j,td) >= 0.d0) then
                    vh(i,j,tc) = v(i,j,tc) + 0.5d0 * phi(i,j,1) * (v(i,j+1,tc) - v(i,j,tc))
                else
                    vh(i,j,tc) = v(i,j+1,tc) + 0.5d0 * phi(i,j,1) * (v(i,j,tc) - v(i,j+1,tc))
                end if
            end do
        end do

        do i = 1,ia+2
            do j = 1,ja+2
                qh(i,j,tc) = (q(i,j,tc) + q(i,j+1,tc)) / 2d0
            end do
        end do

        do i = 1,ia+2
            do j = 1,ja+2
                if (h(i,j,tc) < h_min) then
                    qh(i,j,tc) = 0.d0
                    vh(i,j,tc) = 0.d0
                end if
            end do
        end do
    
    end subroutine cal_vh_and_qh       

    
    subroutine cal_uu(ia,ja,m,n,tc,td,u,q,uu)        ! calculate u and q at control volume u (not the central node)
    
        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) u(ia+3,ja+2,2),q(ia+2,ja+3,2),qud(ia+3,ja+3,2),uu(ia+1,ja+1,2)
        real(kind=8) ruu(m,n+1,1),phi(m,n+1,1),phiruu(ia+1,ja+1,1)

        ! determine the flow direction
        do i = 1,ia+3
            do j = 1,ja+3
                if (i == 1) then
                    qud(i,j,td) = q(i,j,td)
                else if (i == ia + 3) then
                    qud(i,j,td) = q(i-1,j,td)
                else
                    qud(i,j,td) = (q(i-1,j,td) + q(i,j,td)) / 2d0
                end if
            end do
        end do
    
        call cal_rqh(ia,ja,m,n,tc,td,qud,u,ruu)
    
        mm = size(ruu,1)
        nn = size(ruu,2)
        call lim_fun_Q(mm,nn,ruu,phi)
    
        do i = 1,ia+1
            do j = 1,ja+1
                phiruu(i,j,1) = phi(i+1,j+1,1)
                if (qud(i+1,j+1,td) >= 0.d0) then
                    uu(i,j,tc) = u(i+1,j,tc) + 0.5d0 * phiruu(i,j,1) * (u(i+1,j+1,tc) - u(i+1,j,tc))
                else
                    uu(i,j,tc) = u(i+1,j+1,tc) + 0.5d0 * phiruu(i,j,1) * (u(i+1,j,tc) - u(i+1,j+1,tc))
                end if
            end do
        end do

    end subroutine cal_uu
    
    
    subroutine cal_qu(ia,ja,m,n,tc,td,q,uu,qu)

        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) q(ia+2,ja+3,2),uud(ia+3,ja+3,2),uu(ia+1,ja+1,2),qu(ia+1,ja+1,2)
        real(kind=8) rqu(m+1,n,1),phi(m+1,n,1),phirqu(ia+1,ja+1,1)
    
        do i = 1,ia+1
            do j = 1,ja+1
                qu(i,j,tc) = (q(i,j+1,tc) + q(i+1,j+1,tc)) / 2d0
            end do
        end do
    
    end subroutine cal_qu

    
    subroutine cal_vv(ia,ja,m,n,tc,td,v,uu,uud,vv)        ! calculate v and p at control volume v(not the central node)
    
        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) v(ia+2,ja+3,2),uud(ia+3,ja+3,2),vv(ia+1,ja+1,2),uu(ia+1,ja+1,2)
        real(kind=8) rvv(m+1,n,1),phi(m+1,n,1),phirvv(ia+1,ja+1,1)
    
        ! determine the flow direction
        do i = 2,ia+2
            do j = 2,ja+2
                uud(i,j,td) = uu(i-1,j-1,td)         
            end do
        end do
        
        do i = 2,ia+2
            uud(i,1,td) = uu(i-1,1,td)
            uud(i,ja+3,td) = uu(i-1,ja+1,td)
        end do
    
        do j = 1,ja+3
            uud(1,j,td) = uud(2,j,td)
            uud(ia+3,j,td) = uud(ia+2,j,td)
        end do
    
        call cal_rph(ia,ja,m,n,tc,td,uud,v,rvv)
    
        mm = size(rvv,1)
        nn = size(rvv,2)
        call lim_fun_Q(mm,nn,rvv,phi)
    
        do i = 1,ia+1
            do j = 1,ja+1
                phirvv(i,j,1) = phi(i+1,j+1,1)
                if (uud(i+1,j+1,td) >= 0.d0) then
                    vv(i,j,tc) = v(i,j+1,tc) + 0.5d0 * phirvv(i,j,1) * (v(i+1,j+1,tc) - v(i,j+1,tc))
                else
                    vv(i,j,tc) = v(i+1,j+1,tc) + 0.5d0 * phirvv(i,j,1) * (v(i,j+1,tc) - v(i+1,j+1,tc))
                end if
            end do
        end do

    end subroutine cal_vv
    
    
    subroutine cal_pv(ia,ja,m,n,tc,td,p,vv,pv)

        implicit none
    
        integer i,j,ia,ja,m,n,mm,nn,tc,td
        real(kind=8) p(ia+3,ja+2,2),vvd(ia+3,ja+3,2),vv(ia+1,ja+1,2),pv(ia+1,ja+1,2)
        real(kind=8) rpv(m,n+1,1),phi(m,n+1,1),phirpv(ia+1,ja+1,1)
    
        do i = 1,ia+1
            do j = 1,ja+1
                pv(i,j,tc) = (p(i+1,j,1) + p(i+1,j+1,1)) / 2d0
            end do
        end do
    
    end subroutine cal_pv

  
    subroutine cal_deltau(ia,ja,g,ph,uh,qu,uu,hu,eta,deltaweu,deltansu,deltatopou)    ! calculate the x-direction momentum equation 

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) ph(ia+2,ja+2,1),uh(ia+2,ja+2,1),qu(ia+1,ja+1,2),uu(ia+1,ja+1,2)
        real(kind=8) g,hu(ia+3,ja+2,2),eta(ia+2,ja+2,1)
        real(kind=8) deltaweu(ia+1,ja,1),deltansu(ia+1,ja,1),deltatopou(ia+1,ja,1)

        do i = 1,ia+1
            do j = 1,ja
                deltaweu(i,j,1) = ph(i+1,j+1,1) * uh(i+1,j+1,1) - ph(i,j+1,1) * uh(i,j+1,1)
            end do
        end do
    
        do i = 1,ia+1
            do j = 1,ja
                deltansu(i,j,1) = qu(i,j+1,1) * uu(i,j+1,1) - qu(i,j,1) * uu(i,j,1)
            end do
        end do
    
        do i = 1,ia+1
            do j = 1,ja
                deltatopou(i,j,1) = g * hu(i+1,j+1,2) * (eta(i+1,j+1,1) - eta(i,j+1,1))
            end do
        end do
    
    end subroutine cal_deltau
    
    
    subroutine cal_diffu(ia,ja,vis,dx,dy,hu,u,diffu)    

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,vis,hu(ia+3,ja+2,2),u(ia+3,ja+2,2),diffu(ia+1,ja,1)

        do i = 1,ia+1
            do j = 1,ja
                diffu(i,j,1) = vis * (1d0 / (dx)**2d0 * ((hu(i+2,j+1,2) * u(i+2,j+1,1) - &
                    hu(i+1,j+1,2) * u(i+1,j+1,1)) - (hu(i+1,j+1,2) * u(i+1,j+1,1) - &
                    hu(i,j+1,2) * u(i,j+1,1))) + 1d0 / (dy)**2d0 * ((hu(i+1,j+2,2) * u(i+1,j+2,1) - &
                    hu(i+1,j+1,2) * u(i+1,j+1,1)) - (hu(i+1,j+1,2) * u(i+1,j+1,1) - &
                    hu(i+1,j,2) * u(i+1,j,1))))
            end do
        end do
    
    end subroutine cal_diffu

    
    subroutine cal_deltav(ia,ja,g,qh,vh,pv,vv,hv,eta,deltawev,deltansv,deltatopov)    ! calculate the y-direction momentum equation 

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) qh(ia+2,ja+2,1),vh(ia+2,ja+2,1),pv(ia+1,ja+1,2),vv(ia+1,ja+1,2)
        real(kind=8) g,hv(ia+2,ja+3,2),eta(ia+2,ja+2,1)
        real(kind=8) deltawev(ia,ja+1,1),deltansv(ia,ja+1,1),deltatopov(ia,ja+1,1)

        do i = 1,ia
            do j = 1,ja+1
                deltawev(i,j,1) = pv(i+1,j,1) * vv(i+1,j,1) - pv(i,j,1) * vv(i,j,1)
            end do
        end do
    
        do i = 1,ia
            do j = 1,ja+1
                deltansv(i,j,1) = (qh(i+1,j+1,1) * vh(i+1,j+1,1) - qh(i+1,j,1) * vh(i+1,j,1))
            end do
        end do
    
        do i = 1,ia
            do j = 1,ja+1
                deltatopov(i,j,1) = g * hv(i+1,j+1,2) * (eta(i+1,j+1,1) - eta(i+1,j,1))
            end do
        end do
    
    end subroutine cal_deltav    

    
    subroutine cal_diffv(ia,ja,vis,dx,dy,hv,v,diffv)   

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,vis,hv(ia+2,ja+3,2),v(ia+2,ja+3,2),diffv(ia,ja+1,1)

        do i = 1,ia
            do j = 1,ja+1
                diffv(i,j,1) = vis * (1d0 / (dx)**2d0 * ((hv(i+2,j+1,2) * v(i+2,j+1,1) - &
                    hv(i+1,j+1,2) * v(i+1,j+1,1)) - (hv(i+1,j+1,2) * v(i+1,j+1,1) - &
                    hv(i,j+1,2) * v(i,j+1,1))) + 1d0 / (dy)**2d0 * ((hv(i+1,j+2,2) * v(i+1,j+2,1) - &
                    hv(i+1,j+1,2) * v(i+1,j+1,1)) - (hv(i+1,j+1,2) * v(i+1,j+1,1) - &
                    hv(i+1,j,2) * v(i+1,j,1))))
            end do
        end do
    
    end subroutine cal_diffv
    
    
    subroutine cal_u(ia,ja,dx,dy,dt,p_ini,hu,deltaweu,deltansu,deltatopou,diffu,p,u)      ! calculate the new flux in u cv

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,dt,p_ini,hu(ia+3,ja+2,2),u(ia+3,ja+2,2),p(ia+3,ja+2,2)
        real(kind=8) deltaweu(ia+1,ja,1),deltansu(ia+1,ja,1),deltatopou(ia+1,ja,1),diffu(ia+1,ja,1)
    
        do i = 1,ia+1
            do j = 1,ja
                p(i+1,j+1,2) = hu(i+1,j+1,1) * u(i+1,j+1,1) - dt / dx * deltaweu(i,j,1) - &
                    dt / dy * deltansu(i,j,1) - dt / dx * deltatopou(i,j,1) + dt * diffu(i,j,1)
            end do
        end do
    
        ! free flow boundary
        ! set ns boundary conditions      
        do i = 2,ia+2
            p(i,1,2) = p(i,2,2)            
            p(i,ja+2,2) = p(i,ja+1,2)
        end do
        ! set we boundary conditions
        do j = 1,ja+2
            p(1,j,2) = p_ini
            p(ia+3,j,2) = p(ia+2,j,2)
        end do
    
        ! calculate the new velocity in u cv
        do i = 1,ia+3
            do j = 1,ja+2
                if (hu(i,j,2) == 0.d0) then
                    u(i,j,2) = 0.d0
                else
                    u(i,j,2) = p(i,j,2) / hu(i,j,2)
                end if
            end do 
        end do
    
    end subroutine cal_u
 
  
    subroutine cal_v(ia,ja,dx,dy,dt,q_ini,hv,deltawev,deltansv,deltatopov,diffv,q,v)      ! calculate the new flux in v cv

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,dt,q_ini,hv(ia+2,ja+3,2),v(ia+2,ja+3,2),q(ia+2,ja+3,2)
        real(kind=8) deltawev(ia,ja+1,1),deltansv(ia,ja+1,1),deltatopov(ia,ja+1,1),diffv(ia,ja+1,1)
    
        do i = 1,ia
            do j = 1,ja+1
                q(i+1,j+1,2) = hv(i+1,j+1,1) * v(i+1,j+1,1) - dt / dx * deltawev(i,j,1) - &
                    dt / dy * deltansv(i,j,1) - dt / dy * deltatopov(i,j,1) + dt * diffv(i,j,1)
            end do
        end do
    
        ! free flow boundary
        ! set ns boundary conditions
        do i = 2,ia+1
            q(i,1,2) = q(i,2,2)          
            q(i,ja+3,2) = q(i,ja+2,2)
        end do
        ! set we boundary conditions      
        do j = 1,ja+3            
            q(1,j,2) = q_ini
            q(ia+2,j,2) = q(ia+1,j,2)       
        end do
    
        ! calculate the new velocity in u cv
        do i = 1,ia+2
            do j = 1,ja+3
                if (hv(i,j,2) == 0.d0) then
                    v(i,j,2) = 0.d0
                else
                    v(i,j,2) = q(i,j,2) / hv(i,j,2)
                end if
            end do 
        end do
    
    end subroutine cal_v


    subroutine cal_convergence(ia,ja,ia_s,ja_s,u,u_err,u_total)        ! determin convergence

        implicit none
    
        integer i,j,ia,ja,ia_s,ja_s
        real(kind=8) u_err,u_total,u(ia+3,ja+2,2)
    
        u_err = 0.d0
        u_total = 0.d0
    
        do i = 1,ia_s+3
            do j = 1,ja_s+2
                u_err = u_err + abs(u(i,j,2) - u(i,j,1))
                u_total = u_total + abs(u(i,j,1))
            end do
        end do
    
    end subroutine cal_convergence
 
    
    subroutine update_vel(ia,ja,h,u,v)

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) h(ia+2,ja+2,2),u(ia+3,ja+2,2),v(ia+2,ja+3,2)

        ! update the data
        do i = 1,ia+2
            do j = 1,ja+2
                h(i,j,1) = h(i,j,2)
            end do
        end do
    
        do i = 1,ia+3
            do j = 1,ja+2
                u(i,j,1) = u(i,j,2)
            end do
        end do
    
        do i = 1,ia+2
            do j = 1,ja+3
                v(i,j,1) = v(i,j,2)
            end do
        end do
    
    end subroutine update_vel

    
    ! flux limiter parameter r for the water depth in the x-direction    
    subroutine cal_rph(ia,ja,m,n,tc,td,u,h,rph)

        implicit none
    
        integer i,j,ia,ja,m,n,tc,td
        real(kind=8) h(m,n,2),u(m+1,n,2),rph(m+1,n,1)

        do i = 1,m+1
            do j = 1,n
                ! boundary
                if (i == 1 .or. i == m + 1) then       ! ghost boundary
                    rph(i,j,1) = 0.d0
                else if (i == 2) then                ! Leonard mirror node extrapolation
                    if (u(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i-1,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0              
                        else
                            rph(i,j,1) = 2d0 * (h(i-1,j,tc) - h(i-1,j,tc)) / (h(i,j,tc) - h(i-1,j,tc))
                        end if
                    else
                        if (h(i-1,j,tc) - h(i,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0
                        else
                            rph(i,j,1) = (h(i,j,tc) - h(i+1,j,tc)) / (h(i-1,j,tc) - h(i,j,tc))
                        end if
                    end if
                else if (i == m) then
                    if (u(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i-1,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0
                        else
                            rph(i,j,1) = (h(i-1,j,tc) - h(i-2,j,tc)) / (h(i,j,tc) - h(i-1,j,tc))
                        end if
                    else 
                        if (h(i-1,j,tc) - h(i,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0
                        else
                            rph(i,j,1) = 2d0 * (h(i,j,tc) - h(i,j,tc)) / (h(i-1,j,tc) - h(i,j,tc))   
                        end if
                    end if
                else
                    if (u(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i-1,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0
                        else
                            rph(i,j,1) = (h(i-1,j,tc) - h(i-2,j,tc)) / (h(i,j,tc) - h(i-1,j,tc))
                        end if
                    else
                        if (h(i-1,j,tc) - h(i,j,tc) == 0.d0) then
                            rph(i,j,1) = 0.d0
                        else
                            rph(i,j,1) = (h(i,j,tc) - h(i+1,j,tc)) / (h(i-1,j,tc) - h(i,j,tc))
                        end if
                    end if
                end if
            end do
        end do

    end subroutine cal_rph

    
    ! flux limiter parameter r for the velocity in the x-direction
    subroutine cal_rp(ia,ja,m,n,tc,td,phd,u,rp)

        implicit none
    
        integer i,j,ia,ja,tc,td,m,n
        real(kind=8) phd(m-1,n,2),u(m,n,2),rp(m-1,n,1)

        do i = 1,m-1
            do j = 1,n
        ! boundary
                if (i == 1) then                ! Leonard mirror node extrapolation
                    if (phd(i,j,td) >= 0.d0) then
                        if (u(i+1,j,tc) - u(i,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0            
                        else
                            rp(i,j,1) = 2d0 * (u(i,j,tc) - u(i,j,tc)) / (u(i+1,j,tc) - u(i,j,tc))
                        end if
                    else
                        if (u(i,j,tc) - u(i+1,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0
                        else
                            rp(i,j,1) = (u(i+1,j,tc) - u(i+2,j,tc)) / (u(i,j,tc) - u(i+1,j,tc))
                        end if
                    end if
                else if (i == m - 1) then
                    if (phd(i,j,td) >= 0.d0) then
                        if (u(i+1,j,tc) - u(i,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0
                        else
                            rp(i,j,1) = (u(i,j,tc) - u(i-1,j,tc)) / (u(i+1,j,tc) - u(i,j,tc))
                        end if
                    else
                        if (u(i,j,tc) - u(i+1,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0
                        else
                            rp(i,j,1) = 2d0 * (u(i+1,j,tc) - u(i+1,j,tc)) / (u(i,j,tc) - u(i+1,j,tc))    
                        end if
                    end if
                else
                    if (phd(i,j,td) >= 0.d0) then
                        if (u(i+1,j,tc) - u(i,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0
                        else
                            rp(i,j,1) = (u(i,j,tc) - u(i-1,j,tc)) / (u(i+1,j,tc) - u(i,j,tc))
                        end if
                    else
                        if (u(i,j,tc) - u(i+1,j,tc) == 0.d0) then
                            rp(i,j,1) = 0.d0
                        else
                            rp(i,j,1) = (u(i+1,j,tc) - u(i+2,j,tc)) / (u(i,j,tc) - u(i+1,j,tc))
                        end if
                    end if
                end if
            end do
        end do
    
    end subroutine cal_rp
    

    ! flux limiter parameter r for the water depth in the y-direction
    subroutine cal_rqh(ia,ja,m,n,tc,td,v,h,rqh)
    
        implicit none
    
        integer i,j,ia,ja,tc,td,m,n
        real(kind=8) h(m,n,2),v(m,n+1,2),rqh(m,n+1,1)

        do i = 1,m
            do j = 1,n+1
                ! boundary
                if (j == 1 .or. j == n + 1) then          ! ghost boundary
                    rqh(i,j,1) = 0.d0
                else if (j == 2) then                 ! Leonard mirror node extrapolation
                    if (v(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i,j-1,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = 2d0 * (h(i,j-1,tc) - h(i,j-1,tc)) / (h(i,j,tc) - h(i,j-1,tc))
                        end if
                    else
                        if (h(i,j-1,tc) - h(i,j,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = (h(i,j,tc) - h(i,j+1,tc)) / (h(i,j-1,tc) - h(i,j,tc))        
                        end if
                    end if
                else if (j == n) then
                    if (v(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i,j-1,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = (h(i,j-1,tc) - h(i,j-2,tc)) / (h(i,j,tc) - h(i,j-1,tc))
                        end if
                    else
                        if (h(i,j-1,tc) - h(i,j,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = 2d0 * (h(i,j,tc) - h(i,j,tc)) / (h(i,j-1,tc) - h(i,j,tc))       
                        end if
                    end if
                else
                    if (v(i,j,td) >= 0.d0) then
                        if (h(i,j,tc) - h(i,j-1,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = (h(i,j-1,tc) - h(i,j-2,tc)) / (h(i,j,tc) - h(i,j-1,tc))
                        end if
                    else
                        if (h(i,j-1,tc) - h(i,j,tc) == 0.d0) then
                            rqh(i,j,1) = 0.d0
                        else
                            rqh(i,j,1) = (h(i,j,tc) - h(i,j+1,tc)) / (h(i,j-1,tc) - h(i,j,tc))
                        end if
                    end if
                end if
            end do
        end do
    
    end subroutine cal_rqh
    
    
    ! flux limiter parameter r for the velocity in the y-direction
    subroutine cal_rq(ia,ja,m,n,tc,td,qhd,v,rq)

        implicit none
    
        integer i,j,ia,ja,tc,td,m,n
        real(kind=8) qhd(m,n-1,2),v(m,n,2),rq(m,n-1,1)

        do i = 1,m
            do j = 1,n-1
        ! boundary
                if (j == 1) then                ! Leonard mirror node extrapolation
                    if (qhd(i,j,td) >= 0.d0) then
                        if (v(i,j+1,tc) - v(i,j,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0            
                        else
                            rq(i,j,1) = 2d0 * (v(i,j,tc) - v(i,j,tc)) / (v(i,j+1,tc) - v(i,j,tc))
                        end if
                    else
                        if (v(i,j,tc) - v(i,j+1,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0
                        else
                            rq(i,j,1) = (v(i,j+1,tc) - v(i,j+2,tc)) / (v(i,j,tc) - v(i,j+1,tc))
                        end if
                    end if
                else if (j == n - 1) then
                    if (qhd(i,j,td) >= 0.d0) then
                        if (v(i,j+1,tc) - v(i,j,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0
                        else
                            rq(i,j,1) = (v(i,j,tc) - v(i,j-1,tc)) / (v(i,j+1,tc) - v(i,j,tc))
                        end if
                    else
                        if (v(i,j,tc) - v(i,j+1,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0
                        else
                            rq(i,j,1) = 2d0 * (v(i,j+1,tc) - v(i,j+1,tc)) / (v(i,j,tc) - v(i,j+1,tc))   
                        end if
                    end if
                else
                    if (qhd(i,j,td) >= 0.d0) then
                        if (v(i,j+1,tc) - v(i,j,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0
                        else
                            rq(i,j,1) = (v(i,j,tc) - v(i,j-1,tc)) / (v(i,j+1,tc) - v(i,j,tc))
                        end if
                    else
                        if (v(i,j,tc) - v(i,j+1,tc) == 0.d0) then
                            rq(i,j,1) = 0.d0
                        else
                            rq(i,j,1) = (v(i,j+1,tc) - v(i,j+2,tc)) / (v(i,j,tc) - v(i,j+1,tc))
                        end if
                    end if
                end if
            end do
        end do
    
    end subroutine cal_rq
    
    
    ! flux limiter function
    subroutine lim_fun_Q(mm,nn,r,phi)

        implicit none

        integer i,j,mm,nn
        real(kind=8) r(mm,nn,1),phi(mm,nn,1)
        real(kind=8) phi1,phi2,phi3
    
        ! QUICK
        do i = 1,mm
            do j = 1,nn
                phi1 = 2d0 * r(i,j,1)
                phi2 = (3d0 + r(i,j,1)) / 4d0
                phi3 = min(phi1,phi2,2d0)
                phi(i,j,1) = max(0.d0,phi3)
            end do
        end do
    
    end subroutine lim_fun_Q
    
    
    subroutine lim_fun_U(mm,nn,r,phi)

        implicit none

        integer i,j,mm,nn
        real(kind=8) r(mm,nn,1),phi(mm,nn,1)
    
        ! UD
        do i = 1,mm
            do j = 1,nn
                phi(i,j,1) = 0d0
            end do
        end do
    
    end subroutine lim_fun_U
    
    
    subroutine lim_fun_C(mm,nn,r,phi)

        implicit none

        integer i,j,mm,nn
        real(kind=8) r(mm,nn,1),phi(mm,nn,1)
    
        ! CD
        do i = 1,mm
            do j = 1,nn
                phi(i,j,1) = 1d0
            end do
        end do
    
    end subroutine lim_fun_C
    
end module flow_field
