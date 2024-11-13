! topography.f90
!
!****************************************************************************
!
!   This file is part of the Yardang Dynamics Model.
!
!   Set the initial topography.
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


module topography
    
    implicit none
    
    contains
    
    subroutine ini_topography(pi,ia,ja,ia_s,ja_s,lx_s,ly_s,xx_s,yy_s,B)

        implicit none
    
        integer i,j,ia,ja,ia_s,ja_s
        real(kind=8) pi,lx_s,ly_s,h_dune,h_sand,a
        real(kind=8) xx_s(ia_s+2,ja_s+2,1),yy_s(ia_s+2,ja_s+2,1),cir(ia_s+2,ja_s+2,1),B(ia+2,ja+2,2)
    
        h_dune = 2d0      ! dune height
        h_sand = 3d0      ! sand thickness
        a = 50d0          ! dune radius
    
        do i = 1,ia+2
            do j = 1,ja+2
                if (i <= ia_s+2) then
                    ! central coordinate (lx/2,ly/2)
                    cir(i,j,1) = sqrt((xx_s(i,j,1)- lx_s / 2d0)**2d0 + (yy_s(i,j,1) - ly_s / 2d0)**2d0)       
                    if (cir(i,j,1)/a <= 1d0) then
                        B(i,j,1) = h_sand + h_dune * cos(pi * cir(i,j,1) / 2d0 / a)**2d0        ! unit: meters
                    else 
                        B(i,j,1) = h_sand        
                    end if
                else
                    B(i,j,1) = h_sand        
                end if
            end do
        end do
    
    end subroutine ini_topography
    
end module topography