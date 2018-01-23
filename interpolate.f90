module interpolate
    implicit none

    contains
        subroutine interp_linear ( m, data_num, t_data, p_data, interp_num, t_interp, p_interp )
        !*****************************************************************************80
        !
        !! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
        !
        !  Discussion:
        !
        !    From a space of M dimensions, we are given a sequence of
        !    DATA_NUM points, which are presumed to be successive samples
        !    from a curve of points P.
        !
        !    We are also given a parameterization of this data, that is,
        !    an associated sequence of DATA_NUM values of a variable T.
        !    The values of T are assumed to be strictly increasing.
        !
        !    Thus, we have a sequence of values P(T), where T is a scalar,
        !    and each value of P is of dimension M.
        !
        !    We are then given INTERP_NUM values of T, for which values P
        !    are to be produced, by linear interpolation of the data we are given.
        !
        !    Note that the user may request extrapolation.  This occurs whenever
        !    a T_INTERP value is less than the minimum T_DATA or greater than the
        !    maximum T_DATA.  In that case, linear extrapolation is used.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    03 December 2007
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) M, the spatial dimension.
        !
        !    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
        !
        !    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
        !    independent variable at the sample points.  The values of T_DATA
        !    must be strictly increasing.
        !
        !    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
        !    dependent variables at the sample points.
        !
        !    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
        !    at which interpolation is to be done.
        !
        !    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
        !    independent variable at the interpolation points.
        !
        !    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
        !    values of the dependent variables at the interpolation points.
        !

        integer ( kind = 4 ) data_num
        integer ( kind = 4 ) m
        integer ( kind = 4 ) interp_num

        integer ( kind = 4 ) interp
        integer ( kind = 4 ) left
        real ( kind = 8 ) p_data(m,data_num)
        real ( kind = 8 ) p_interp(m,interp_num)
        logical r8vec_ascends_strictly
        integer ( kind = 4 ) right
        real ( kind = 8 ) t
        real ( kind = 8 ) t_data(data_num)
        real ( kind = 8 ) t_interp(interp_num)

          do interp = 1, interp_num

            t = t_interp(interp)
        !
        !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
        !  nearest to, TVAL.
        !
            call r8vec_bracket ( data_num, t_data, t, left, right )

            p_interp(1:m,interp) = ((t_data(right)-t)*p_data(1:m,left)+(t &
                - t_data(left))*p_data(1:m,right))/(t_data(right)-t_data(left))

          end do

          return
        end subroutine

        subroutine r8vec_bracket ( n, x, xval, left, right )
        !*****************************************************************************80
        !
        !! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
        !
        !  Discussion:
        !
        !    An R8VEC is an array of double precision real values.
        !
        !    If the values in the vector are thought of as defining intervals
        !    on the real line, then this routine searches for the interval
        !    nearest to or containing the given value.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    06 April 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, length of input array.
        !
        !    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
        !
        !    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
        !
        !    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
        !    Either:
        !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
        !      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
        !    or
        !      X(LEFT) <= XVAL <= X(RIGHT).
        !
        integer ( kind = 4 ) n
        integer ( kind = 4 ) i
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) xval
        integer ( kind = 4 ) left
        integer ( kind = 4 ) right

          do i = 2, n - 1

            if ( xval < x(i) ) then
              left = i - 1
              right = i
              return
            end if

           end do

          left = n - 1
          right = n

          return
        end subroutine
end module interpolate