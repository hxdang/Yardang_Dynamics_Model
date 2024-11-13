!  Yardang_Dynamics_Model.f90 
!
!****************************************************************************
!
!   This is the main program for the Yardang Dynamics Model.
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


program Yardang_Dynamics_Model

    use grid
    use topography
    use ini_condition
    use flow_field
    use sand_transport

    implicit none

    include 'parameter_flow.txt'      ! read calculation parameters
    include 'parameter_sand.txt'
    
    integer i,j,tc,td,m,n
    integer num,selectedGrids,ri,rj,ti,tj,ava_num,dep_est,count_num,seed(1)
    real(kind=8) t,dt,tstart,tend_vel,tend_sur,tend,t_vel,t_sur,t_total,dtmin,dt_sur,C,u_err,u_total
    real(kind=8) dx,dy,xx,yy,B,h,eta,deltat,nu,h_in,h_out
    real(kind=8) u,uu,uh,p,pv,Qsx,hu,huc,ph,uud,puv,uv
    real(kind=8) v,vv,vh,q,qu,Qsy,hv,hvc,qh,qvu,vu
    real(kind=8) deltaweu,deltansu,deltatopou,diffu,deltawev,deltansv,deltatopov,diffv
    real(kind=8) uf,vf,vel,velf,tan_thetax,tan_thetay,tan_theta,theta_vel
    real(kind=8) BR,S,VEG,slope,uft,vft,velft,qsat,slope_value,ero
    real(kind=8) hu_s(ia_s+3,ja_s+2,1),u_s(ia_s+3,ja_s+2,1),hv_s(ia_s+2,ja_s+3,1),v_s(ia_s+2,ja_s+3,1)
    real(kind=8) xx_s(ia_s+2,ja_s+2,1),yy_s(ia_s+2,ja_s+2,1),thetarep(ia_s+2,ja_s+2,1)
    character(20) fileNumberStr,fileNumber
    
    common /time /t,dt,dt_sur
    common /CFL /deltat(ia,ja,1),nu(ia,ja,1)
    common /mesh /dx,dy,xx(ia+2,ja+2,1),yy(ia+2,ja+2,1)
    common /solution /B(ia+2,ja+2,2),h(ia+2,ja+2,2),eta(ia+2,ja+2,1),hu(ia+3,ja+2,2),hv(ia+2,ja+3,2),&
                      u(ia+3,ja+2,2),p(ia+3,ja+2,2),Qsx(ia+1,ja,1),v(ia+2,ja+3,2),q(ia+2,ja+3,2),Qsy(ia,ja+1,1),&
                      ph(ia+2,ja+2,1),uh(ia+2,ja+2,1),qh(ia+2,ja+2,1),vh(ia+2,ja+2,1),uud(ia+3,ja+3,2),&
                      qu(ia+1,ja+1,2),uu(ia+1,ja+1,2),pv(ia+1,ja+1,2),vv(ia+1,ja+1,2),h_in(ja+2,1),h_out(ja+2,1),&
                      huc(ia+3,ja+2,2),deltaweu(ia+1,ja,1),deltansu(ia+1,ja,1),deltatopou(ia+1,ja,1),diffu(ia+1,ja,1),&
                      hvc(ia+2,ja+3,2),deltawev(ia,ja+1,1),deltansv(ia,ja+1,1),deltatopov(ia,ja+1,1),diffv(ia,ja+1,1),&
                      qvu(ia+1,ja,1),vu(ia+1,ja,1),puv(ia,ja+1,1),uv(ia,ja+1,1),selectedGrids(ia_s,ja_s),&
                      uf(ia_s+3,ja_s+2,1),vf(ia_s+2,ja_s+3,1),velf(ia_s+2,ja_s+2,1),&
                      tan_thetax(ia_s+1,ja_s,2),tan_thetay(ia_s,ja_s+1,2),vel(ia_s+2,ja_s+2,1),theta_vel(ia+2,ja_s+2,1),&
                      uft(ia_s+1,ja_s,2),vft(ia_s,ja_s+1,2),qsat(ia_s+2,ja_s+2,1),tan_theta(ia_s+2,ja_s+2,1),&
                      BR(ia_s+2,ja_s+2,2),S(ia_s+2,ja_s+2,1),VEG(ia_s+2,ja_s+2,1),slope(ia_s+2,ja_s+2,2)
    
    dx = lx / ia
    dy = ly / ja
    
    tstart = 0.d0
    
    dtmin = 0.03d0
    tend_vel = 300d0
    
    dt_sur = 60d0 * 60d0 * 2d0              ! time of an iteration
    tend_sur = 60d0 * 60d0 * 24d0 * 10d0          ! duration of an evolution
    tend = 60d0 * 60d0 * 24d0 * 10d0        ! total running time
    
    
    call griding(ia,ja,dx,dy,xx,yy)
    
    call griding(ia_s,ja_s,dx,dy,xx_s,yy_s)
    
    call ini_topography(pi,ia,ja,ia_s,ja_s,lx_s,ly_s,xx_s,yy_s,B)
    
    ! 将初始地形数据写入文件
    open(unit = 10,file = 'result\ini_topography.txt')
    write(10,*) 'xx               yy                 B'
    do i = 1,ia_s+2
        do j = 1,ja_s+2
            write(10,*) xx_s(i,j,1),yy_s(i,j,1),B(i,j,1)
        end do
    end do
    close(10)
    
    t_total = tstart
    num = 1

    call random_seed()
    
    do while(t_total <= tend)
        
    call set_zero(ia,ja,h,u,v,p,q,qu,uu,pv,vv,hu,hv)

    call initial(h_ini,h_min,p_ini,q_ini,ia,ja,B,h,p,q,hu,hv,u,v)
    
    t_vel = tstart
    u_err = 1d0
    u_total = 1d0
    
    ! flow calculation start
    do while (t_vel <= tend_vel)
    
    m = size(h,1)
    n = size(h,2)
    tc = 1          ! select the time step for the calculation
    td = 1
    call cal_hu(ia,ja,m,n,tc,td,h_min,h,u,hu)
    call cal_hv(ia,ja,m,n,tc,td,h_min,h,v,hv)
    call cal_huc(ia,ja,tc,h_min,h,huc)
    call cal_hvc(ia,ja,tc,h_min,h,hvc)
    
    call cal_p(ia,ja,p_ini,u,hu,p)
    
    call cal_q(ia,ja,q_ini,v,hv,q)
    
    call cal_dt_and_C(ia,ja,dx,dy,g,dtmin,t_vel,u,v,p,q,h,dt,C)
    
    call cal_h(ia,ja,h_min,dt,dx,dy,h_in,h_out,h,p,q)
    
    call cal_eta(ia,ja,h,B,eta)
    
    tc = 2
    td = 1
    call cal_hu(ia,ja,m,n,tc,td,h_min,h,u,hu)
    call cal_hv(ia,ja,m,n,tc,td,h_min,h,v,hv)
    call cal_huc(ia,ja,tc,h_min,h,huc)
    call cal_hvc(ia,ja,tc,h_min,h,hvc)
    
    m = size(u,1)
    n = size(u,2)
    tc = 1
    td = 1
    call cal_uh_and_ph(ia,ja,m,n,tc,td,h_min,h,u,p,uh,ph)
    
    m = size(v,1)
    n = size(v,2)
    tc = 1
    td = 1
    call cal_vh_and_qh(ia,ja,m,n,tc,td,h_min,h,v,q,vh,qh)
    
    m = size(u,1)
    n = size(u,2)
    tc = 1
    td = 1
    call cal_uu(ia,ja,m,n,tc,td,u,q,uu)
    
    m = size(q,1)
    n = size(q,2)
    tc = 1
    td = 1
    call cal_qu(ia,ja,m,n,tc,td,q,uu,qu)
    
    m = size(v,1)
    n = size(v,2)
    tc = 1
    td = 1
    call cal_vv(ia,ja,m,n,tc,td,v,uu,uud,vv)
    
    m = size(p,1)
    n = size(p,2)
    tc = 1
    td = 1
    call cal_pv(ia,ja,m,n,tc,td,p,vv,pv)
    
    call cal_deltau(ia,ja,g,ph,uh,qu,uu,huc,eta,deltaweu,deltansu,deltatopou)
    
    call cal_diffu(ia,ja,vis,dx,dy,huc,u,diffu)
    
    call cal_u(ia,ja,dx,dy,dt,p_ini,huc,deltaweu,deltansu,deltatopou,diffu,p,u)
    
    call cal_deltav(ia,ja,g,qh,vh,pv,vv,hvc,eta,deltawev,deltansv,deltatopov)
    
    call cal_diffv(ia,ja,vis,dx,dy,hvc,v,diffv)
    
    call cal_v(ia,ja,dx,dy,dt,q_ini,hvc,deltawev,deltansv,deltatopov,diffv,q,v)
         
    call cal_convergence(ia,ja,ia_s,ja_s,u,u_err,u_total)
    
    if (u(ia+2,(ja+2)/2,2) > u(ia,(ja+2)/2,2)) then
        exit
    end if
    
    call update_vel(ia,ja,h,u,v)
    
    t_vel = t_vel + dt
    
    end do
    
    !print *,'The computation reaches convergence when t = ',t_vel
    print *,'reach convergence when t =',t_vel
    print *,'relative error =',u_err/u_total
    
    call update_vel(ia,ja,h,u,v)
    
    call transmit_data(ia,ja,ia_s,ja_s,hu,hv,u,v,hu_s,hv_s,u_s,v_s)   ! transmit data
    
    call initialize(ia,ja,ia_s,ja_s,lx_s,ly_s,pi,xx_s,yy_s,B,thetarep,BR,S,VEG)
    
    call cal_vel(ia_s,ja_s,u_s,v_s,theta_vel,vel)
    
    call cal_frivel(ia_s,ja_s,karman,z0,hu_s,hv_s,u_s,v_s,uf,vf,velf)
    
    call cal_angle_of_slope(ia_s,ja_s,dx,dy,S,tan_thetax,tan_thetay,tan_theta)
    
    call cal_critical_frivel(ia_s,ja_s,fricoe,g,ds,rho_air,s_rho,thetarep,tan_thetax,tan_thetay,uft,vft,velft)
    
    call cal_qsat(ia_s,ja_s,g,rho_air,velf,velft,qsat)
    
    ! evolution start
    t_sur = tstart
    
    do while (t_sur < tend_sur)
        
        ava_num = n_total
        selectedGrids = 0
        do while (ava_num > 0)
            
            call choose_erosion_grid(ia_s,ja_s,ava_num,selectedGrids,ri,rj)
            ava_num = ava_num - 1
            
            call cal_erosion(ia_s,ja_s,ri,rj,dx,dy,dt_sur,rho_sand,thetarep,VEG,qsat,theta_vel,ero,slope,S)
            
            dep_est = 0
            ti = ri
            tj = rj
            do while (dep_est == 0)
                
                if (ero == 0.d0) then
                    exit
                end if
                
                call cal_transport(ia_s,ja_s,dx,dy,g,vis,rho_sand,rho_air,ds,velft,theta_vel,velf,u_s,v_s,ti,tj)
                
                call cal_deposition(ia_s,ja_s,ti,tj,dx,dy,thetarep,ero,velft,velf,slope,VEG,dep_est,S)
                
            end do
            
        end do
            
        t_sur = t_sur + dt_sur
    end do
    
    call update_topo(ia,ja,ia_s,ja_s,S,B)
    
    t_total = t_total + t_sur
    
    ! save the topography data
    write(fileNumberStr, '(I0)') int(t_total/3600d0)
    open(unit = 1,file = 'result\topography_t=' // trim(adjustl(fileNumberStr)) // 'h.txt')
    write(1,*) 'xx               yy                 B'
    do i = 1,ia_s+2
        do j = 1,ja_s+2
            write(1,*) xx_s(i,j,1),yy_s(i,j,1),B(i,j,1)
        end do
    end do
    close(1)
    
    ! save the velocity data
    open(unit = 2,file = 'result\u_t=' // trim(adjustl(fileNumberStr)) // 'h.txt')
    write(2,*) 'hu               u'
    do i = 1,ia_s+3
        do j = 1,ja_s+2
            write(2,*) hu(i,j,2),u(i,j,2)
        end do
    end do
    close(2)
    
    open(unit = 3,file = 'result\v_t=' // trim(adjustl(fileNumberStr)) // 'h.txt')
    write(3,*) 'hv               v'
    do i = 1,ia_s+2
        do j = 1,ja_s+3
            write(3,*) hv(i,j,2),v(i,j,2)
        end do
    end do
    close(3)
    
    print *,num,'times evolution finished'
    num = num + 1
    
    end do
    
    print *,"Finish!"
    
end program Yardang_Dynamics_Model