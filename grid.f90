! grid.f90
!
!****************************************************************************
!
!   This file is part of the Yardang Dynamics Model.
!
!   Generate the grid.
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


module grid
    
    implicit none
    
    contains
    
    subroutine griding(ia,ja,dx,dy,xx,yy)

        implicit none
    
        integer i,j,ia,ja
        real(kind=8) dx,dy,xx(ia+2,ja+2,1),yy(ia+2,ja+2,1)
    
        do i = 1,ia+2
            do j = 1,ja+2
                xx(i,j,1) = -dx + (i - 0.5d0) * dx
                yy(i,j,1) = -dy + (j - 0.5d0) * dy
            end do
        end do
    
    end subroutine griding
    
end module grid