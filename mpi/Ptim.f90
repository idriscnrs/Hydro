!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Ptim.f90 --- Timing MPI or mixed MPI/OpenMP Fortran programs
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Thu Nov 18 15:44:01 2004
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Wed Aug  4 16:32:28 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE PTIM
  IMPLICIT NONE

!  PRIVATE

  !... Shared private variables
  INTEGER, PARAMETER         :: ROOT_RANK=0, d=SELECTED_REAL_KIND(12)
  REAL(kind=d)               :: ELA_overhead, CPU_overhead
  REAL(kind=d), DIMENSION(2) :: ELA_tim, CPU_tim
  LOGICAL                    :: mpi=.FALSE.

  !$ INTEGER, PARAMETER                           :: MAX_NUM_THRDS=128
  !$ REAL(kind=d), DIMENSION(0:MAX_NUM_THRDS-1,2) :: ELA_omp_tim=0.0_d
  !$ INTEGER                                      :: num_thrds=1
  !$ LOGICAL                                      :: omp=.FALSE.

  INTEGER, DIMENSION(8)           :: values
  CHARACTER(LEN=8), DIMENSION(2)  :: date
  CHARACTER(LEN=10), DIMENSION(2) :: time
  CHARACTER(LEN=5)                :: zone
  CHARACTER(LEN=80)               :: user_label="My Program"
  INTEGER                         :: label_length=10

  !... Public interfaces
  PUBLIC :: PTIM_start, PTIM_stop
  !$ PUBLIC :: PTIM_omp_start, PTIM_omp_stop

  CONTAINS

  SUBROUTINE PTIM_start(label)
    IMPLICIT NONE

    !... MPI Header file
    INCLUDE "mpif.h"

    !... Input dummy argument
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: label

    !... Local variables
    INTEGER      :: rank, ierr
    REAL(kind=d) :: dummy

    !... Compute clock overhead
    ELA_overhead = MPI_WTIME()
    ELA_overhead = MPI_WTIME() - ELA_overhead
    CALL CPU_TIME(dummy)
    CALL CPU_TIME(CPU_overhead)
    IF (dummy < 0.0_d) PRINT *,"Warning, PTIM_start: CPU user time is not available on this machine."
    CPU_overhead = CPU_overhead - dummy

    Call MPI_INITIALIZED(mpi, ierr)
    if (.not. mpi) PRINT *,"Warning, PTIM_start: MPI_init must be invoked prior to PTIM_start."
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    !... Start timings on "date & time"
    IF ( rank == ROOT_RANK ) THEN
       CALL DATE_AND_TIME(date(1),time(1),zone,values)
       IF(PRESENT(label)) THEN
         label_length=LEN_TRIM(label)
         user_label(1:label_length) = TRIM(label)
       END IF
    END IF
    !... Start elapsed and CPU time counters
    ELA_tim(1) = MPI_WTIME()
    CALL CPU_TIME(CPU_tim(1))
  END SUBROUTINE PTIM_start

  !$ SUBROUTINE PTIM_omp_start()
    !$ USE OMP_LIB
    !$ IMPLICIT NONE

    !... Local variables
    !$ INTEGER :: rank, nt

    !$ omp=OMP_IN_PARALLEL()
    !$ IF (.NOT. omp) PRINT *,"Warning, PTIM_omp_start must be called in an OpenMP PARALLEL scope."

    !$ nt = OMP_GET_NUM_THREADS()
    !$ num_thrds=MAX(num_thrds, nt)
    !$ IF (num_thrds > MAX_NUM_THRDS) THEN
    !$   PRINT *,"Warning, PTIM_omp_start: Number of threads greater than MAX_NUM_THDS=128."
    !$   PRINT *,"Please, reduce the number of threads in the PARALLEL scope or set MAX_NUM_THDS"
    !$   PRINT *,"to a greater value in PTIM module."
    !$ END IF

    ! Start thread time counter.
    !$ rank = OMP_GET_THREAD_NUM()
    !$ ELA_omp_tim(rank,1)=OMP_GET_WTIME()
  !$ END SUBROUTINE PTIM_omp_start

  !$ SUBROUTINE PTIM_omp_stop()
    !$ USE OMP_LIB
    !$ IMPLICIT NONE

    !... Local variables
    !$ INTEGER      :: rank
    !$ REAL(KIND=d) :: omp_fin_elapse

    ! Store final private time counter.
    !$ omp_fin_elapse=OMP_GET_WTIME()

    !$ omp=OMP_IN_PARALLEL()
    !$ IF (.NOT. omp) PRINT *,"Warning, PTIM_omp_stop must be called in an OpenMP PARALLEL scope."

    ! Which thread I am ?
    !$ rank=OMP_GET_THREAD_NUM()

    ! Collecting each private elapsed time into a global array.
    !$ ELA_omp_tim(rank,2) = ELA_omp_tim(rank,2) + omp_fin_elapse - ELA_omp_tim(rank,1)
  !$ END SUBROUTINE PTIM_omp_stop

  SUBROUTINE PTIM_stop(label, file)
    IMPLICIT NONE

    !... MPI Header file
    INCLUDE "mpif.h"

    !... Input dummy arguments
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: label, file

    !... Local variables
    INTEGER                                 :: rank, nb_procs, ierr, length, i, u=6
    LOGICAL                                 :: file_opened
    REAL(kind=d), ALLOCATABLE, DIMENSION(:) :: all_ELA_tim, all_CPU_tim, all_ratio
    REAL(kind=d)                            :: tot_ELA_tim, tot_CPU_tim, tot_ratio,&
                                               max_ELA_tim, max_CPU_tim, max_ratio,&
                                               min_ELA_tim, min_CPU_tim, min_ratio,&
                                               avg_ELA_tim, avg_CPU_tim, avg_ratio

    !$ INTEGER, ALLOCATABLE, DIMENSION(:)        :: all_num_thrds
    !$ LOGICAL, ALLOCATABLE, DIMENSION(:)        :: all_omp
    !$ REAL(kind=d), ALLOCATABLE, DIMENSION(:,:) :: all_ELA_omp_tim
    !$ INTEGER                                   :: j, proc_max_num_thrds
    

    CHARACTER(LEN=256) :: prologue,epilogue,proc_tab_label,proc_tab_2hline,&
                          proc_tab_fmt
    !$ CHARACTER(LEN=256) :: omp_tab_label,omp_tab_hline,omp_tab_fmt, omp_tab_longhline
    CHARACTER(LEN=512)  :: others

    !... Final CPU and elapsed times
    CALL CPU_TIME(CPU_tim(2))
    ELA_tim(2) = MPI_WTIME() - ELA_tim(1) - ELA_overhead - CPU_overhead
    CPU_tim(2) = CPU_tim(2) - CPU_tim(1) - CPU_overhead

    !... Gather all timings
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, ierr)
    IF ( rank == ROOT_RANK ) THEN
      IF (PRESENT(label)) THEN
        length=LEN_TRIM(label)
        IF (VERIFY(label(1:length), user_label(1:label_length)) /= 0) &
           PRINT *,"Warning, PTIM_stop: Input label argument doesn't match the one of PTIM_start argument."
      END IF
      IF (PRESENT(file)) THEN
        !... Looking for a none connected logical file unit.
        u=99
        INQUIRE(UNIT=u, OPENED=file_opened)
        DO WHILE (file_opened .AND. u /= 0)
          u = u - 1
          INQUIRE(UNIT=u, OPENED=file_opened)
        END DO
        IF(u == 0 .AND. file_opened) THEN
          PRINT *,"Warning, PTIM_stop: All file units from 0 to 99 are already connected. Output has been redirected to STDOUT."
          u=6
        ELSE
          OPEN(UNIT=u, FILE=TRIM(file), FORM="formatted", STATUS="unknown", ACTION="write")
        END IF
      END IF        

      ALLOCATE(all_ELA_tim(0:nb_procs-1), all_CPU_tim(0:nb_procs-1), all_ratio(0:nb_procs-1))
      !$ ALLOCATE(all_omp(0:nb_procs-1))
    END IF

    CALL MPI_GATHER(ELA_tim(2), 1, MPI_DOUBLE_PRECISION,  &
                    all_ELA_tim, 1, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(CPU_tim(2), 1, MPI_DOUBLE_PRECISION,  &
                    all_CPU_tim, 1, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    !$ CALL MPI_GATHER(omp, 1, MPI_LOGICAL, all_omp, 1, MPI_LOGICAL, &
    !$               & 0, MPI_COMM_WORLD, ierr)

    !$ ALLOCATE(all_num_thrds(0:nb_procs-1))
    !$ CALL MPI_ALLGATHER(num_thrds, 1, MPI_INTEGER, all_num_thrds, 1, MPI_INTEGER, &
    !$                  & MPI_COMM_WORLD, ierr)
    !$ proc_max_num_thrds=MAXVAL(all_num_thrds(0:))

    !$ IF (rank == ROOT_RANK ) ALLOCATE(all_ELA_omp_tim(0:proc_max_num_thrds-1,0:nb_procs-1))
    !$ CALL MPI_GATHER(ELA_omp_tim(0,2), proc_max_num_thrds, MPI_DOUBLE_PRECISION, &
    !$               & all_ELA_omp_tim, proc_max_num_thrds, MPI_DOUBLE_PRECISION,  &
    !$               & 0, MPI_COMM_WORLD, ierr)

    IF ( rank == ROOT_RANK) THEN
      !... Compute elapsed user time
       tot_ELA_tim = SUM(all_ELA_tim(0:))
       avg_ELA_tim = tot_ELA_tim/REAL(nb_procs,kind=d)
       max_ELA_tim = MAXVAL(all_ELA_tim(0:))
       min_ELA_tim = MINVAL(all_ELA_tim(0:))
       IF ( min_ELA_tim <= 0.0_d ) THEN
          PRINT *,"Warning, PTIM_stop: Measured elapsed user time seems to be too short"
          PRINT *,"compared to the clock precision. Timings could be erroneous."
       END IF

       !... Compute CPU user time
       tot_CPU_tim = SUM(all_CPU_tim(0:))
       avg_CPU_tim = tot_CPU_tim/REAL(nb_procs,kind=d)
       max_CPU_tim = MAXVAL(all_CPU_tim(0:))
       min_CPU_tim = MINVAL(all_CPU_tim(0:))
       IF ( min_CPU_tim <= 0.0_d ) THEN
          PRINT *,"Warning, PTIM_stop: Measured CPU user time seems to be too short"
          PRINT *,"compared to the clock precision. Timings could be erroneous."
       END IF

       !... Compute cpu/elapsed ratio
       all_ratio(0:) = all_CPU_tim(0:) / all_ELA_tim(0:)
       tot_ratio     = SUM(all_ratio(0:))
       min_ratio     = MINVAL(all_ratio(0:))
       max_ratio     = MAXVAL(all_ratio(0:))
       avg_ratio     = tot_ratio/REAL(nb_procs,kind=d)

       !... End of timings at "date & time"
       CALL DATE_AND_TIME(date(2),time(2),zone,values)

       !... Output Format
       prologue       ='(//,5X,"Copyright (C) 2004, Jalel CHERGUI, IDRIS-CNRS FRANCE.",/,&
                                         &5X,"TIM (release 3.1) summary report of *** ",A," ***",/)'
       proc_tab_label ='(5X,"Process Rank",(" "),"|",(" "),"Process CPU Time (s)",(" "),"|",&
                                      &(" "),"Process Elapsed Time (s)",(" "),"|",(" "),"CPU/Elapsed")'
       !$ omp_tab_label   ='(5X,13(" "),"|",22(" "),"|",(" "),"Thread Elapsed Time (s)",2(" "),"|")'
       proc_tab_2hline='(5X,13("="),"|",22("="),"|",26("="),"|",12("="))'
       !$ omp_tab_longhline ='(5X,13("-"),"|",22("-"),"|",26("-"),"|",12("-"))'
       !$ omp_tab_hline   ='(5X,13(" "),"|",22(" "),"|",26("-"),"|")'
       proc_tab_fmt   ='(5X,I4,9(" "),"|",3(" "),F12.3,7(" "),"|",4(" "),F12.3,10(" "),"|",5(" "),F7.3)'
       !$ omp_tab_fmt     ='(5X,13(" "),"|",22(" "),"|",(" "),I4,":",4(" "),F12.3,4(" "),"|")'
       others         ='(5X,"Total",  8(" "),"|",3(" "),F12.3,7(" "),"|",4(" "),F12.3,10(" "),"|",5(" "),F7.3,/,&
                                      &5X,"Minimum",6(" "),"|",3(" "),F12.3,7(" "),"|",4(" "),F12.3,10(" "),"|",5(" "),F7.3,/,&
                                      &5X,"Maximum",6(" "),"|",3(" "),F12.3,7(" "),"|",4(" "),F12.3,10(" "),"|",5(" "),F7.3,/,&
                                      &5X,"Average",6(" "),"|",3(" "),F12.3,7(" "),"|",4(" "),F12.3,10(" "),"|",5(" "),F7.3)'
       epilogue       ='(/,5X,"TIM started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT",/,&
                                      &5X,  "TIM stopped on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT",/)'
       WRITE(UNIT=u,FMT=TRIM(prologue)) user_label(1:label_length)
       WRITE(UNIT=u,FMT=TRIM(proc_tab_label))
       WRITE(UNIT=u,FMT=TRIM(proc_tab_2hline))
       ! DO i = 0, nb_procs-1
       !   !$ IF (all_omp(i) .AND. all_num_thrds(i) > 1) THEN
       !   !$   WRITE(UNIT=u,FMT=TRIM(proc_tab_fmt)) i, all_CPU_tim(i), &
       !   !$                                        all_ELA_tim(i), all_ratio(i)
       !   !$   WRITE(UNIT=u,FMT=TRIM(omp_tab_hline))
       !   !$   WRITE(UNIT=u,FMT=TRIM(omp_tab_label))
       !   !$   WRITE(UNIT=u,FMT=TRIM(omp_tab_hline))
       !   !$   DO j = 0, all_num_thrds(i)-1
       !   !$     WRITE(UNIT=u,FMT=TRIM(omp_tab_fmt)) j, all_ELA_omp_tim(j,i)
       !   !$   END DO
       !   !$     IF ( i < nb_procs-1) WRITE(UNIT=u,FMT=TRIM(omp_tab_longhline))
       !   !$ ELSE IF (.not.all_omp(i)) THEN
       !     WRITE(UNIT=u,FMT=TRIM(proc_tab_fmt)) i, all_CPU_tim(i), &
       !                                          all_ELA_tim(i), all_ratio(i)
       !   !$ END IF
       ! END DO
       WRITE(UNIT=u,FMT=TRIM(proc_tab_2hline))
       WRITE(UNIT=u,FMT=TRIM(others)) tot_CPU_tim, tot_ELA_tim, tot_ratio, &
                                      min_CPU_tim, min_ELA_tim, min_ratio, &
                                      max_CPU_tim, max_ELA_tim, max_ratio, &
                                      avg_CPU_tim, avg_ELA_tim, avg_ratio
       WRITE(UNIT=u,FMT=TRIM(proc_tab_2hline))
       WRITE(UNIT=u,FMT=TRIM(epilogue)) date(1)(7:8), date(1)(5:6), date(1)(1:4), &
                                        time(1)(1:2), time(1)(3:4), time(1)(5:6), &
                                        zone(1:3),    zone(4:5),                  &
                                        date(2)(7:8), date(2)(5:6), date(2)(1:4), &
                                        time(2)(1:2), time(2)(3:4), time(2)(5:6), &
                                        zone(1:3),    zone(4:5)
       IF(PRESENT(file) .AND. u /= 6) CLOSE(UNIT=u)
       !$ DEALLOCATE(all_ELA_omp_tim, all_omp)
       DEALLOCATE(all_ELA_tim, all_CPU_tim, all_ratio)
     END IF
     !$ DEALLOCATE(all_num_thrds)

     !... RAZ
     !$ ELA_omp_tim(:,:)=0.0_d ; num_thrds=1 ; omp=.false.
  END SUBROUTINE PTIM_stop
END MODULE PTIM
