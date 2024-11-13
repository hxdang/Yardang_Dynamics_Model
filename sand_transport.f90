! sand_transport.f90
!
!****************************************************************************
!
!   This file is part of the Yardang Dynamics Model.
!
!   Calculte the sand transport.
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

    
module sand_transport
    
    implicit none
    
    contains
    
    subroutine transmit_data(ia,ja,ia_s,ja_s,hu,hv,u,v,hu_s,hv_s,u_s,v_s)          

        implicit none
    
        integer ia,ja,ia_s,ja_s
        real(kind=8) hu(ia+3,ja+2,2),u(ia+3,ja+2,2),hv(ia+2,ja+3,2),v(ia+2,ja+3,2)
        real(kind=8) hu_s(ia_s+3,ja_s+2,1),u_s(ia_s+3,ja_s+2,1),hv_s(ia_s+2,ja_s+3,1),v_s(ia_s+2,ja_s+3,1)
    
        hu_s(1:ia_s+3,1:ja_s+2,1) = hu(1:ia_s+3,1:ja_s+2,2)
        hv_s(1:ia_s+2,1:ja_s+3,1) = hv(1:ia_s+2,1:ja_s+3,2)
        u_s(1:ia_s+3,1:ja_s+2,1) = u(1:ia_s+3,1:ja_s+2,2)
        v_s(1:ia_s+2,1:ja_s+3,1) = v(1:ia_s+2,1:ja_s+3,2)
    
    end subroutine transmit_data   
    
    
    subroutine initialize(ia,ja,ia_s,ja_s,lx_s,ly_s,pi,xx_s,yy_s,B,thetarep,BR,S,VEG)            

        implicit none
    
        integer i,j,ia,ja,ia_s,ja_s
        real(kind=8) lx_s,ly_s,pi
        real(kind=8) xx_s(ia_s+2,ja_s+2,1),yy_s(ia_s+2,ja_s+2,1),B(ia+2,ja+2,2),cir(ia_s+2,ja_s+2,1)
        real(kind=8) BR(ia_s+2,ja_s+2,2),S(ia_s+2,ja_s+2,1),VEG(ia_s+2,ja_s+2,1),thetarep(ia_s+2,ja_s+2,1)
    
        do i = 1,ia_s+2
            do j = 1,ja_s+2
                S(i,j,1) = B(i,j,1)
                VEG(i,j,1) = 0.d0       ! vegetation cover
                thetarep(i,j,1) = 32d0 * pi / 180d0     ! angle of repose
            end do
        end do
    
    end subroutine initialize   
  
    
    subroutine cal_vel(ia,ja,u,v,theta_vel,vel)         

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) u(ia+3,ja+2,1),v(ia+2,ja+3,1),vel(ia+2,ja+2,1),tantheta_vel(ia+2,ja+2,1),theta_vel(ia+2,ja+2,1)
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (u(i+1,j,1) == 0.d0) then
                    tantheta_vel(i,j,1) = huge(tantheta_vel(i,j,1))
                else
                    tantheta_vel(i,j,1) = abs(v(i,j+1,1)) / abs(u(i+1,j,1))              ! angle between velocity and x-axis
                end if
                theta_vel(i,j,1) = atan(tantheta_vel(i,j,1))
                vel(i,j,1) = sqrt(u(i+1,j,1)**2d0 + v(i,j+1,1)**2d0)        
            end do
        end do
    
    end subroutine cal_vel
    
    
    subroutine cal_frivel(ia,ja,karman,z0,hu,hv,u,v,uf,vf,velf)            ! calculate the friction velocity

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) z0,karman,hu(ia+3,ja+2,1),u(ia+3,ja+2,1),hv(ia+2,ja+3,1),v(ia+2,ja+3,1)
        real(kind=8) uf(ia+3,ja+2,1),vf(ia+2,ja+3,1),velf(ia+2,ja+2,1)
    
        do i = 1,ia+3
            do j = 1,ja+2
                uf(i,j,1) = karman * u(i,j,1) / log(hu(i,j,1) / z0)
                if (u(i,j,1) > 0.d0 .and. u(i,j,1) < 1d-6) then           
                    uf(i,j,1) = 0.d0
                elseif (u(i,j,1) < 0.d0 .and. u(i,j,1) > -1d-6) then
                    uf(i,j,1) = 0.d0
                end if
            end do
        end do
    
        do i = 1,ia+2
            do j = 1,ja+3
                vf(i,j,1) = karman * v(i,j,1) / log(hv(i,j,1) / z0)
                if (v(i,j,1) > 0.d0 .and. v(i,j,1) < 1d-6) then
                    vf(i,j,1) = 0.d0
                elseif (v(i,j,1) < 0.d0 .and. v(i,j,1) > -1d-6) then
                    vf(i,j,1) = 0.d0
                end if
            end do
        end do
    
        do i = 1,ia+2
            do j = 1,ja+2
                velf(i,j,1) = sqrt(uf(i+1,j,1)**2d0 + vf(i,j+1,1)**2d0)
            end do
        end do
    
    end subroutine cal_frivel
    
    
    subroutine cal_angle_of_slope(ia,ja,dx,dy,S,tan_thetax,tan_thetay,tan_theta)            ! calculate the angle of slope

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,S(ia+2,ja+2,1),tan_thetax(ia+1,ja,2),tan_thetay(ia,ja+1,2),tan_theta(ia+2,ja+2,1)
    
        do i = 1,ia+1
            do j = 1,ja
                tan_thetax(i,j,1) = (S(i+1,j+1,1) - S(i,j+1,1)) / dx
            end do
        end do
    
        do i = 1,ia
            do j = 1,ja+1
                tan_thetay(i,j,1) = (S(i+1,j+1,1) - S(i+1,j,1)) / dy
            end do
        end do
    
        do i = 2,ia+1
            do j = 2,ja+1
                tan_theta(i,j,1) = sqrt(((S(i+1,j,1) - S(i-1,j,1)) / 2d0 / dx)**2d0 + ((S(i,j+1,1) - S(i,j-1,1)) / 2d0 / dy)**2d0)
            end do
        end do
    
        ! set ns boundary conditions
        do i = 2,ia+1
            tan_theta(i,1,1) = tan_theta(i,2,1)           
            tan_theta(i,ja+2,1) = tan_theta(i,ja+1,1)
        end do
        ! set we boundary conditions
        do j = 1,ja+2       
            tan_theta(1,j,1) = tan_theta(2,j,1)
            tan_theta(ia+2,j,1) = tan_theta(ia+1,j,1)
        end do
    
    end subroutine cal_angle_of_slope  
    

    ! calculate the critical friction velocity
    subroutine cal_critical_frivel(ia,ja,fricoe,g,ds,rho_air,s_rho,thetarep,tan_thetax,tan_thetay,uft,vft,velft) 

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) fricoe,g,ds,rho_air,rho_sand,s_rho,thetarep(ia+2,ja+2,1),A_t,velft
        real(kind=8) thetax,sin_thetax,cos_thetax,thetay,sin_thetay,cos_thetay
        real(kind=8) tan_thetax(ia+1,ja,2),tan_thetay(ia,ja+1,2),uft(ia+1,ja,2),vft(ia,ja+1,2)
    
        do i = 1,ia+1
            do j = 1,ja
                thetax = atan(abs(tan_thetax(i,j,1)))
                sin_thetax = sin(thetax)
                cos_thetax = cos(thetax)
                uft(i,j,1) = fricoe * sqrt((rho_sand - rho_air) * g * ds * &
                            (cos_thetax + sin_thetax / tan(thetarep(i,j+1,1))) / rho_air)
            end do
        end do
    
        do i = 1,ia
            do j = 1,ja+1
                thetay = atan(abs(tan_thetay(i,j,1)))
                sin_thetay = sin(thetay)
                cos_thetay = cos(thetay)
                vft(i,j,1) = fricoe * sqrt((rho_sand - rho_air) * g * ds * &
                            (cos_thetay + sin_thetay / tan(thetarep(i+1,j,1))) / rho_air)
            end do
        end do
    
        A_t = 0.1d0
        velft = A_t * sqrt(ds * g * (s_rho - 1d0))
    
    end subroutine cal_critical_frivel
    
    
    subroutine cal_qsat(ia,ja,g,rho_air,velf,velft,qsat) 

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) C_sat,g,rho_air,velft,velf(ia+2,ja+2,1),qsat(ia+2,ja+2,1)
    
        C_sat = 2.78d0
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (abs(velf(i,j,1)) <= velft) then
                    qsat(i,j,1) = 0.d0
                else
                    qsat(i,j,1) = C_sat * rho_air / g * (velf(i,j,1) - velft) * (velf(i,j,1) + velft)**2d0
                end if
            end do
        end do
    
    end subroutine cal_qsat
    
    
    subroutine choose_erosion_grid(ia,ja,ava_num,selectedGrids,ri,rj)           ! select an erosion grid randomly

        implicit none
    
        integer size,i,j,ia,ja,n_total,count,ava_num,availableGrids(ava_num)
        integer randomIndex,min_value,max_value,selectedGrids(ia,ja)
        integer col,row,ri,rj
    
        ! list of selectable grids
        count = 0
        do i = 1,ia
            do j = 1,ja
                if (selectedGrids(i,j) == 0) then
                  count = count + 1
                  availableGrids(count) = j + ja * (i - 1)
                end if
            end do
        end do 

        ! select a grid randomly 
        min_value = 1
        max_value = ava_num
        call choose_range_random_number(min_value,max_value,randomIndex)
    
        row = mod(availableGrids(randomIndex)-1,ja) + 1
        col = (availableGrids(randomIndex) - row) / ja + 1
        rj = row + 1
        ri = col + 1
    
        selectedGrids(col,row) = 1
    
    end subroutine choose_erosion_grid

    
    subroutine cal_erosion(ia,ja,ri,rj,dx,dy,dt_sur,rho_sand,thetarep,VEG,qsat,theta_vel,ero,slope,S)           ! erosion process

        implicit none
    
        integer i,j,ia,ja,ri,rj,steep_i,steep_in,steep_j,steep_jn
        real(kind=8) dx,dy,dt_sur,rho_sand,thetarep(ia+2,ja+2,1),ero,p_ero,slope_value,random_num
        real(kind=8) S(ia+2,ja+2,1),slope(ia+2,ja+2,2),VEG(ia+2,ja+2,1),qsat(ia+2,ja+2,1),theta_vel(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)
        
        if (S(ri,rj,1) < 1d-6) then
            S(ri,rj,1) = S(ri,rj,1)
            ero = 0.d0
        else
            p_ero = 1d0 - VEG(ri,rj,1)

            call choose_random_number(random_num)

            if (random_num > p_ero) then
                S(ri,rj,1) = S(ri,rj,1)
                ero = 0.d0
            else
                ero = qsat(ri,rj,1) * dt_sur * (dx / cos(theta_vel(ri,rj,1))) / rho_sand / (dx * dy)
                S(ri,rj,1) = S(ri,rj,1) - ero
                ! slope limitation
                call cal_slope(ri,rj,ia,ja,dx,dy,S,slope_value)
                slope(ri,rj,1) = slope_value
    
                if (slope(ri,rj,1) <= tan(thetarep(ri,rj,1))) then
                    S(ri,rj,1) = S(ri,rj,1)
                else
                    call cal_steepest_ero(ri,rj,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                    call choose_steepest_grid(ri,rj,local_height,top_three(1),steep_i,steep_j)
                
                    if (steep_j == 0) then
                        S(ri,rj,1) = S(ri,rj,1)
                        steep_jn = 0
                    else if (steep_j == -1) then
                        S(ri,rj,1) = S(ri,rj,1)
                        steep_jn = 0
                    else
                        call slope_correction(ri,rj,steep_i,steep_j,ia,ja,dx,dy,thetarep,S)
                    
                        call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                        slope(steep_i,steep_j,1) = slope_value
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    end if
                
                    do while (steep_jn > 0)
                    
                        call slope_correction(steep_i,steep_j,steep_in,steep_jn,ia,ja,dx,dy,thetarep,S)
                    
                        ! calculate the steepest grid of the new target grid 
                        steep_i = steep_in
                        steep_j = steep_jn

                        call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                        slope(steep_i,steep_j,1) = slope_value
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    
                    end do
                end if
            end if
        end if
    
    end subroutine cal_erosion
    
    
    subroutine cal_transport(ia,ja,dx,dy,g,vis,rho_sand,rho_air,ds,velft,theta_vel,velf,u,v,ti,tj) 

        implicit none
    
        integer i,j,ia,ja,di,dj,ti,tj
        real(kind=8) dx,dy,g,vis,l_a,l_b,rho_sand,ds,rho_air,velft,l_sat,l_satx,l_saty
        real(kind=8) u(ia+3,ja+2,1),v(ia+2,ja+3,1),theta_vel(ia+2,ja+2,1),velf(ia+2,ja+2,1)
    
        l_sat = 8
        l_satx =  l_sat * cos(theta_vel(ti,tj,1))
        l_saty =  l_sat * sin(theta_vel(ti,tj,1))
    
        ! transport direction
        if (u(ti,tj,1) >= 0.d0) then
            di = 0
        else
            di = 1
        end if
        if (v(ti,tj,1) >= 0.d0) then
            dj = 0
        else
            dj = 1
        end if
    
        ! target grid 
        ti = ti + (-1)**di * int(l_satx / dx + 0.5d0)         
        tj = tj + (-1)**dj * int(l_saty / dy + 0.5d0)
    
        ! periodic boundary
        if (ti > ia + 1) then
            ti = 1 + ti - (ia + 1)
        elseif (ti < 2) then
            ti = ia + ti
        end if
        if (tj > ja + 1) then
            tj = 1 + tj - (ja + 1)
        elseif (tj < 2) then
            tj = ja + tj
        end if
    
    end subroutine cal_transport
    
    
    subroutine cal_deposition(ia,ja,ti,tj,dx,dy,thetarep,ero,velft,velf,slope,VEG,dep_est,S)           ! deposition process

        implicit none
    
        integer i,j,ia,ja,ti,tj,dep_est,rep_i1,rep_j1,rep_i2,rep_j2,steep_i,steep_j,steep_in,steep_jn
        real(kind=8) dx,dy,thetarep(ia+2,ja+2,1),velft
        real(kind=8) ero,rep,p_sp,p_dep,p_rep,slope_value,delta,rep1,rep2,random_num
        real(kind=8) velf(ia+2,ja+2,1),S(ia+2,ja+2,1),slope(ia+2,ja+2,2),VEG(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)
    
        ! probability of deposition 
        if (S(ti,tj,1) < 1d-6) then
            p_sp = 0.4d0
        else
            p_sp = 0.6d0
        end if

        if (velf(ti,tj,1) < velft) then
            p_dep = 1d0
        else
            p_dep = p_sp + VEG(ti,tj,1) * (1d0 - p_sp)
        end if

        ! probability of reptation
        p_rep = 1 - VEG(ti,tj,1)

        rep = ero

        call choose_random_number(random_num)
    
        if (random_num < p_dep) then
            dep_est = 1                           ! deposit
            S(ti,tj,1) = S(ti,tj,1) + ero
        
            ! reptation
            call choose_random_number(random_num)
            if (random_num > p_rep) then        
                S(ti,tj,1) = S(ti,tj,1)       ! no reptation
                delta = ero
            else
                call cal_steepest_reptation(ti,tj,ia,ja,dx,dy,S,local_height,top_three)
                call choose_reptation_grid(ti,tj,local_height,top_three(1),rep_i1,rep_j1)
                call choose_reptation_grid(ti,tj,local_height,top_three(2),rep_i2,rep_j2)
            
                if (rep_j1 == 0 .and. rep_j2 == 0) then                      
                    S(ti,tj,1) = S(ti,tj,1)
                    delta = ero
            
                elseif (rep_j1 > 0 .and. rep_j2 == 0) then
                    S(ti,tj,1) = S(ti,tj,1) - rep
                    S(rep_i1,rep_j1,1) = S(rep_i1,rep_j1,1) + rep
                    delta = ero - rep              
                
                    call slope_limit(rep_i1,rep_j1,ia,ja,dx,dy,thetarep,rep,slope,S,steep_i,steep_j)
                
                else
                    rep1 = rep * top_three(1) / (top_three(1) + top_three(2))
                    rep2 = rep * top_three(2) / (top_three(1) + top_three(2))
                
                    S(ti,tj,1) = S(ti,tj,1) - rep
                    S(rep_i1,rep_j1,1) = S(rep_i1,rep_j1,1) + rep1
                    S(rep_i2,rep_j2,1) = S(rep_i2,rep_j2,1) + rep2
                    delta = ero - rep
                
                    call slope_limit(rep_i1,rep_j1,ia,ja,dx,dy,thetarep,rep1,slope,S,steep_i,steep_j)
                
                    call slope_limit(rep_i2,rep_j2,ia,ja,dx,dy,thetarep,rep2,slope,S,steep_i,steep_j)
                
                end if
            end if
    
            ! slope limitation for target grid N
            call cal_slope(ti,tj,ia,ja,dx,dy,S,slope_value)
            slope(ti,tj,1) = slope_value
        
            if (slope(ti,tj,1) <= tan(thetarep(ti,tj,1))) then
                S(ti,tj,1) = S(ti,tj,1)
            else
                ! erosion or deposition
                if (delta > 0) then        ! deposition
                    call cal_steepest_dep(ti,tj,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                    call choose_steepest_grid(ti,tj,local_height,top_three(1),steep_i,steep_j)
                else
                    call cal_steepest_ero(ti,tj,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                    call choose_steepest_grid(ti,tj,local_height,top_three(1),steep_i,steep_j)
                end if
            
                if (steep_j == 0) then                    
                    S(ti,tj,1) = S(ti,tj,1)
                    steep_jn = 0
                else if (steep_j == -1) then  
                    S(ti,tj,1) = S(ti,tj,1)
                    steep_jn = 0
                else
                    call slope_correction(ti,tj,steep_i,steep_j,ia,ja,dx,dy,thetarep,S)
                    
                    call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                    slope(steep_i,steep_j,1) = slope_value
                    
                    if (delta > 0) then        ! deposition
                        call cal_steepest_dep(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    else
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    end if    
                end if
            
                do while (steep_jn > 0)
                
                    call slope_correction(steep_i,steep_j,steep_in,steep_jn,ia,ja,dx,dy,thetarep,S)
                
                    steep_j = steep_jn
                    steep_i = steep_in
                
                    call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                    slope(steep_i,steep_j,1) = slope_value
                
                    if (delta > 0) then           ! deposition
                        call cal_steepest_dep(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    else
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    end if
                end do  
            end if
        else
            dep_est = 0                     ! bounce
            S(ti,tj,1) = S(ti,tj,1)
        
            ! reptation
            call choose_random_number(random_num)
        
            if (random_num > p_rep) then
                S(ti,tj,1) = S(ti,tj,1)       ! no reptation
                delta = 0
            else
                call cal_steepest_reptation(ti,tj,ia,ja,dx,dy,S,local_height,top_three)
                call choose_reptation_grid(ti,tj,local_height,top_three(1),rep_i1,rep_j1)
                call choose_reptation_grid(ti,tj,local_height,top_three(2),rep_i2,rep_j2)
            
                if (rep_j1 == 0 .and. rep_j2 == 0) then                   
                    S(ti,tj,1) = S(ti,tj,1)
                    delta = 0
                
                elseif (rep_j1 > 0 .and. rep_j2 == 0) then
                    S(ti,tj,1) = S(ti,tj,1) - rep
                    S(rep_i1,rep_j1,1) = S(rep_i1,rep_j1,1) + rep
                    delta = rep
                
                    call slope_limit(rep_i1,rep_j1,ia,ja,dx,dy,thetarep,rep,slope,S,steep_i,steep_j)
                
                else
                    rep1 = rep * top_three(1) / (top_three(1) + top_three(2))
                    rep2 = rep * top_three(2) / (top_three(1) + top_three(2))
                
                    S(ti,tj,1) = S(ti,tj,1) - rep
                    S(rep_i1,rep_j1,1) = S(rep_i1,rep_j1,1) + rep1
                    S(rep_i2,rep_j2,1) = S(rep_i2,rep_j2,1) + rep2
                    delta = rep
                
                    call slope_limit(rep_i1,rep_j1,ia,ja,dx,dy,thetarep,rep1,slope,S,steep_i,steep_j)
                
                    call slope_limit(rep_i2,rep_j2,ia,ja,dx,dy,thetarep,rep2,slope,S,steep_i,steep_j)
                
                end if
            end if
            ! slope limitation for target grid N
            if (delta > 0) then
                call cal_slope(ti,tj,ia,ja,dx,dy,S,slope_value)
                slope(ti,tj,1) = slope_value
            
                if (slope(ti,tj,1) <= tan(thetarep(ti,tj,1))) then
                    S(ti,tj,1) = S(ti,tj,1)
                else
                    call cal_steepest_ero(ti,tj,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                    call choose_steepest_grid(ti,tj,local_height,top_three(1),steep_i,steep_j)
                
                    if (steep_j == 0) then                 
                        S(ti,tj,1) = S(ti,tj,1)
                        steep_jn = 0
                    else if (steep_j == -1) then  
                        S(ti,tj,1) = S(ti,tj,1)
                        steep_jn = 0
                    else
                        call slope_correction(ti,tj,steep_i,steep_j,ia,ja,dx,dy,thetarep,S)
                    
                        call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                        slope(steep_i,steep_j,1) = slope_value
                    
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                    end if
                
                    do while (steep_jn > 0)
                    
                        call slope_correction(steep_i,steep_j,steep_in,steep_jn,ia,ja,dx,dy,thetarep,S)
                    
                        steep_j = steep_jn
                        steep_i = steep_in
                    
                        call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                        slope(steep_i,steep_j,1) = slope_value
                    
                        call cal_steepest_ero(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                        call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)  
                    
                    end do
                end if
            end if
        end if   
    
    end subroutine cal_deposition


    subroutine update_topo(ia,ja,ia_s,ja_s,S,B)

        implicit none
    
        integer i,j,ia,ja,ia_s,ja_s
        real(kind=8) S(ia_s+2,ja_s+2,1),B(ia+2,ja+2,2)

        ! update the data
    
        do i = 1,ia_s+2
            do j = 1,ja_s+2
                B(i,j,1) = S(i,j,1)
            end do
        end do
    
    end subroutine update_topo
    
    
    subroutine cal_slope(i_sl,j_sl,ia,ja,dx,dy,S,slope_value)          ! calculate the slope

        implicit none
    
        integer i,j,i_sl,j_sl,ia,ja
        real(kind=8) dx,dy,slpoe_x,slpoe_y,slope_value,max_dh,local_height(3,3,1),S(ia+2,ja+2,1)
    
        do i = -1,1
            do j = -1,1
                if (j_sl+j == 0 .or. i_sl+i == 0 .or. j_sl+j == ja+3 .or. i_sl+i == ia+3) then
                    local_height(i+2,j+2,1) = -1d5
                else if (i == 0 .and. j == 0) then
                    local_height(i+2,j+2,1) = -1d5
                else
            	    local_height(i+2,j+2,1) = abs(S(i_sl,j_sl,1) - S(i_sl+i,j_sl+j,1))
                end if
            end do
        end do
   
        max_dh = maxval(local_height)
        slpoe_x = max_dh / dx
        slpoe_y = max_dh / dy
        
        slope_value = sqrt(slpoe_x**2d0 + slpoe_y**2d0)

    end subroutine cal_slope
    
    
    subroutine slope_correction(ri,rj,steep_i,steep_j,ia,ja,dx,dy,thetarep,S)         

        implicit none
    
        integer i,j,ia,ja,ri,rj,steep_i,steep_j
        real(kind=8) dx,dy,dh,delta_dh,slope_value,thetarep(ia+2,ja+2,1)
        real(kind=8) local_height(3,3,1),S(ia+2,ja+2,1),slope(ia+2,ja+2,2)
        
            if (S(steep_i,steep_j,1) - S(ri,rj,1) > 0) then
                dh = sqrt(tan(thetarep(ri,rj,1))**2d0 / (1 / dx**2d0 + 1 / dy**2d0))
                delta_dh = (S(steep_i,steep_j,1) - S(ri,rj,1) - dh) / 2d0
                S(ri,rj,1) = S(ri,rj,1) + delta_dh
                S(steep_i,steep_j,1) = S(steep_i,steep_j,1) - delta_dh
            else
                dh = sqrt(tan(thetarep(ri,rj,1))**2d0 / (1 / dx**2d0 + 1 / dy**2d0))
                delta_dh = (S(ri,rj,1) - S(steep_i,steep_j,1) - dh) / 2d0
                S(ri,rj,1) = S(ri,rj,1) - delta_dh
                S(steep_i,steep_j,1) = S(steep_i,steep_j,1) + delta_dh
            end if

    end subroutine slope_correction
    
    
    ! determin the steepest direction for erosion
    subroutine cal_steepest_ero(i_st,j_st,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)          

        implicit none
    
        integer i,j,i_st,j_st,ia,ja
        real(kind=8) dx,dy,thetarep(ia+2,ja+2,1),slope(ia+2,ja+2,2),steep,S(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)
    
        if (slope(i_st,j_st,1) > tan(thetarep(i_st,j_st,1))) then
            do i = -1,1
                do j = -1,1
                    if (j_st+j == 0 .or. i_st+i == 0 .or. j_st+j == ja+3 .or. i_st+i == ia+3) then
                        local_height(i+2,j+2,1) = -1d5
                    elseif (i == 0 .and. j == 0) then
                        local_height(i+2,j+2,1) = -1d5
                    else
                        if (S(i_st+i,j_st+j,1) > S(i_st,j_st,1)) then
                            local_height(i+2,j+2,1) = S(i_st+i,j_st+j,1) - S(i_st,j_st,1)
                        else
                            local_height(i+2,j+2,1) = -1d5
                        end if
                    end if
                end do
            end do
        
            top_three = [local_height(1,1,1),local_height(1,1,1),local_height(1,1,1)]
        
            do i = 1,3
                do j = 1,3
                    if (local_height(i,j,1) > top_three(1)) then
                        top_three(3) = top_three(2)
                        top_three(2) = top_three(1)
                        top_three(1) = local_height(i,j,1)
                    elseif (local_height(i,j,1) > top_three(2)) then
                        top_three(3) = top_three(2)
                        top_three(2) = local_height(i,j,1)
                    elseif (local_height(i,j,1) > top_three(3)) then
                        top_three(3) = local_height(i,j,1)
                    end if
                end do
            end do
        else
            top_three = [1d6,1d6,1d6]
        end if

    end subroutine cal_steepest_ero
    
    
    ! determin the steepest direction for reptation
    subroutine cal_steepest_reptation(i_st2,j_st2,ia,ja,dx,dy,S,local_height,top_three)         

        implicit none
    
        integer i,j,i_st2,j_st2,ia,ja
        real(kind=8) dx,dy,S(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)
    
        do i = -1,1
            do j = -1,1
                if (j_st2+j == 0 .or. i_st2+i == 0 .or. j_st2+j == ja+3 .or. i_st2+i == ia+3) then
                    local_height(i+2,j+2,1) = -1d5
                elseif (i == 0 .and. j == 0) then
                    local_height(i+2,j+2,1) = -1d5
                else
                    if (S(i_st2+i,j_st2+j,1) < S(i_st2,j_st2,1)) then
                        local_height(i+2,j+2,1) = S(i_st2,j_st2,1) - S(i_st2+i,j_st2+j,1)
                    else
                        local_height(i+2,j+2,1) = -1d5
                    end if
                end if
            end do
        end do
        
        top_three = [local_height(1,1,1),local_height(1,1,1),local_height(1,1,1)]

        do i = 1,3
            do j = 1,3
                if (local_height(i,j,1) > top_three(1)) then
                    top_three(3) = top_three(2)
                    top_three(2) = top_three(1)
                    top_three(1) = local_height(i,j,1)
                elseif (local_height(i,j,1) > top_three(2)) then
                    top_three(3) = top_three(2)
                    top_three(2) = local_height(i,j,1)
                elseif (local_height(i,j,1) > top_three(3)) then
                    top_three(3) = local_height(i,j,1)
                end if
            end do
        end do
    
    end subroutine cal_steepest_reptation

    
    ! determin the steepest direction for deposition
    subroutine cal_steepest_dep(i_st,j_st,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)          

        implicit none
    
        integer i,j,i_st,j_st,ia,ja
        real(kind=8) dx,dy,thetarep(ia+2,ja+2,1),slope(ia+2,ja+2,2),S(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)
    
        if (slope(i_st,j_st,1) > tan(thetarep(i_st,j_st,1))) then
            do i = -1,1
                do j = -1,1
                    if (j_st+j == 0 .or. i_st+i == 0 .or. j_st+j == ja+3 .or. i_st+i == ia+3) then
                        local_height(i+2,j+2,1) = -1d5
                    elseif (i == 0 .and. j == 0) then
                        local_height(i+2,j+2,1) = -1d5
                    else
                        if (S(i_st+i,j_st+j,1) < S(i_st,j_st,1)) then
                            local_height(i+2,j+2,1) = S(i_st,j_st,1) - S(i_st+i,j_st+j,1)
                        else
                            local_height(i+2,j+2,1) = -1d5
                        end if
                    end if
                end do
            end do
        
            top_three = [local_height(1,1,1),local_height(1,1,1),local_height(1,1,1)]
            
            do i = 1,3
                do j = 1,3
                    if (local_height(i,j,1) > top_three(1)) then
                        top_three(3) = top_three(2)
                        top_three(2) = top_three(1)
                        top_three(1) = local_height(i,j,1)
                    elseif (local_height(i,j,1) > top_three(2)) then
                        top_three(3) = top_three(2)
                        top_three(2) = local_height(i,j,1)
                    elseif (local_height(i,j,1) > top_three(3)) then
                        top_three(3) = local_height(i,j,1)
                    end if
                end do
            end do
        else
            top_three = [1d6,1d6,1d6]
        end if
    
    end subroutine cal_steepest_dep


    subroutine slope_limit(rep_i,rep_j,ia,ja,dx,dy,thetarep,rep,slope,S,steep_i,steep_j)
    
        implicit none
    
        integer rep_i,rep_j,ia,ja,steep_i,steep_j,steep_in,steep_jn
        real(kind=8) dx,dy,rep,slope_value,thetarep(ia+2,ja+2,1),slope(ia+2,ja+2,2),S(ia+2,ja+2,1)
        real(kind=8) top_three(3),local_height(3,3,1)

        call cal_slope(rep_i,rep_j,ia,ja,dx,dy,S,slope_value)
        slope(rep_i,rep_j,1) = slope_value
                
        if (slope(rep_i,rep_j,1) <= tan(thetarep(rep_i,rep_j,1))) then
            S(rep_i,rep_j,1) = S(rep_i,rep_j,1)
        else
            call cal_steepest_dep(rep_i,rep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
            call choose_steepest_grid(rep_i,rep_j,local_height,top_three(1),steep_i,steep_j)
                    
            if (steep_j == 0) then                     
                S(rep_i,rep_j,1) = S(rep_i,rep_j,1)
                steep_jn = 0
            else if (steep_j == -1) then  
                S(rep_i,rep_j,1) = S(rep_i,rep_j,1)
                steep_jn = 0
            else
            
                call slope_correction(rep_i,rep_j,steep_i,steep_j,ia,ja,dx,dy,thetarep,S)
                
                call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                slope(steep_i,steep_j,1) = slope_value
                call cal_steepest_dep(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)
                        
            end if
                    
            do while (steep_jn > 0)
            
                call slope_correction(steep_i,steep_j,steep_in,steep_jn,ia,ja,dx,dy,thetarep,S)                    
                    
                ! calculate the steepest grid of the new target grid 
                steep_i = steep_in;
                steep_j = steep_jn;
                    
                call cal_slope(steep_i,steep_j,ia,ja,dx,dy,S,slope_value)
                slope(steep_i,steep_j,1) = slope_value
                
                call cal_steepest_dep(steep_i,steep_j,ia,ja,dx,dy,thetarep,slope,S,local_height,top_three)
                call choose_steepest_grid(steep_i,steep_j,local_height,top_three(1),steep_in,steep_jn)

            end do  
        end if
                        
    end subroutine slope_limit

    
    subroutine choose_steepest_grid(i_st,j_st,local_height,steep_value,steep_i,steep_j)

        implicit none
    
        integer i,j,i_st,j_st,col,row,steep_i,steep_j
        real(kind=8) steep_value,local_height(3,3,1)

        if (steep_value == 1d6) then
            steep_j = 0
            steep_i = 0
        else if (steep_value == -1d5) then
            steep_j = -1
            steep_i = -1
        else
            do i = 1,3
                do j = 1,3
                    if (local_height(i,j,1) == steep_value) then
                        col = i
                        row = j
                    end if
                end do
            end do
    
            steep_j = j_st + row - 2
            steep_i = i_st + col - 2
        end if

    end subroutine choose_steepest_grid
    
    
    subroutine choose_reptation_grid(i_st,j_st,local_height,steep_value,steep_i,steep_j)

        implicit none
    
        integer i,j,i_st,j_st,col,row,steep_i,steep_j
        real(kind=8) steep_value,local_height(3,3,1)

        if (steep_value < 0.d0) then
            steep_j = 0
            steep_i = 0
        else
            do i = 1,3
                do j = 1,3
                    if (local_height(i,j,1) == steep_value) then
                        col = i
                        row = j
                    end if
                end do
            end do
    
            steep_j = j_st + row - 2
            steep_i = i_st + col - 2
        end if

    end subroutine choose_reptation_grid
    
    
    subroutine choose_range_random_number(min_value,max_value,randomIndex)

        implicit none
    
        integer randomIndex,min_value,max_value
        real(kind=8) random_value

        call random_number(random_value)
    
        randomIndex = min_value + int(real(max_value - min_value + 1) * random_value)

    end subroutine choose_range_random_number
    
    
    subroutine choose_random_number(random_num)

        implicit none
    
        real(kind=8) random_num

        call random_number(random_num)

    end subroutine choose_random_number
    
end module sand_transport