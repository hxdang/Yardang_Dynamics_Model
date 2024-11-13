! ini_condition.f90
!
!****************************************************************************
!
!   This file is part of the Yardang Dynamics Model.
!
!   Set initial conditions.
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

    
module ini_condition
    
    implicit none
    
    contains
    
    subroutine initial(h_ini,h_min,p_ini,q_ini,ia,ja,B,h,p,q,hu,hv,u,v)

        integer i,ia,j,ja
        real(kind=8) h_ini,h_min,p_ini,q_ini,B(ia+2,ja+2,2),h(ia+2,ja+2,2),h_in(ja+2,1),h_out(ja+2,1)
        real(kind=8) p(ia+3,ja+2,2),hu(ia+3,ja+2,2),u(ia+3,ja+2,2),q(ia+2,ja+3,2),hv(ia+2,ja+3,2),v(ia+2,ja+3,2)
    
        ! initial water depth
        do i = 1,ia+2
            do j = 1,ja+2
                h(i,j,1) = h_ini - B(i,j,1)
            end do
        end do
        
        do i = 1,ia+2
            do j = 1,ja+2
                if (h(i,j,1) < h_min) then
                    h(i,j,1) = 0.d0
                end if
            end do
        end do
    
        do j = 1,ja+2
            h_in(j,1) = h(1,j,1)
            h_out(j,1) = h(ia+2,j,1)
        end do

        ! initial flux
        do i = 1,1
            do j = 1,ja+2
                p(i,j,1) = p_ini
            end do
        end do
    
        do i = 1,1
            do j = 1,ja+3
                q(i,j,1) = q_ini
            end do
        end do

        ! initial velocity
        do i = 1,1
            do j = 1,ja+2
                u(i,j,1) = p_ini / h_in(j,1)
            end do
        end do
    
        do i = 1,1
            do j = 1,ja+3
                if (j == ja+3) then
                    v(i,j,1) = q_ini / h_in(j-1,1)
                else
                    v(i,j,1) = q_ini / h_in(j,1)
                end if
            end do
        end do
    
    end subroutine initial
    
end module ini_condition