       use eesunhong_astropy_interface,
     &     only : outer_multiply
       use stdlib_kinds, only : sp
       EXTERNAL FCN
       real(sp) :: outer_product
       outer_product = outer_multiply(4.5, 3.0)
       CALL MINUIT(FCN,0)
       END
