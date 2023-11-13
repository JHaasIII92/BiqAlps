      subroutine dtrsl(t, ldt, n, b, job, info)
      integer ldt, n, job, info
      double precision t(ldt,1), b(1)
       character*1 uplo, trans
 
       if (job .eq. 00) then
           uplo = 'L'
           trans = 'N'
       else if (job .eq. 01) then
           uplo = 'U'
           trans = 'N'
       else if (job .eq. 10) then
           uplo = 'L'
           trans = 'T'
       else if (job .eq. 11) then
           uplo = 'U'
           trans = 'T'
       endif
       call dtrtrs(uplo, trans, 'N', n, 1, t, ldt, b, n, info)
       end
c====================== The end of dtrsl ==============================
      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,1)
      call dpotrf('U', n, a, lda, info)
      end
c====================== The end of dpofa ===============================
