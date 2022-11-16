c
c
c
c
      module cl_option
      implicit none 
      character*40 argv, cl_outfile
      integer*4 idx
      integer*8 cl_seed, cl_numEvents, cl_type
      
      contains 
        subroutine read_commandPara(numOpts)
        integer*4 numOpts
        idx = 1
        do while (idx < numOpts)
          call get_command_argument(idx, argv)
          select case (argv)
             case ("--seed")
                 idx = idx + 1
                 call get_command_argument(idx, argv)
                 read(argv, *) cl_seed
                 write(*,*) "seed number: ", cl_seed
             case ("--output")
                 idx = idx + 1
                 call get_command_argument(idx, argv)
                 cl_outfile = argv
                 write(*,*) "output file: ", cl_outfile
             case ("--nevents")
                 idx = idx + 1
                 call get_command_argument(idx, argv)
                 read(argv, *) cl_numEvents
                 write(*,*) "number of events: ", cl_numEvents
             case ("--type")
                 idx = idx + 1
                 call get_command_argument(idx, argv)
                 read(argv, *) cl_type
                 if (cl_type .eq. 1) then
                    write(*,*) "type: ", cl_type, 
     &                         "  (radiated inclusive inelastics)"
                 else if (cl_type .eq. 2) then
                    write(*,*) "type: ", cl_type, 
     &                         "  (radiated elastics)"
                 else if (cl_type .eq. 3) then
                    write(*,*) "type: ", cl_type, 
     &                         "  (radiated QE)"
                 else
                    write(*,*) "Wrong type!"
                 endif
             case default 
                 write(*,*) "Invaild option used, ignore."
          end select
          idx = idx + 1
c           write(*,*) idx, argv
        enddo
    
        end subroutine read_commandPara   

      end module cl_option
