       use eesunhong_astropy_interface,
     &     only : outer_multiply
       EXTERNAL FCN
       call outer_multiply(4.5, 3.0)
       CALL MINUIT(FCN,0)
       END
