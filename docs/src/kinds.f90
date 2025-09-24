module hyperSIS_kinds_mod
    use, intrinsic :: iso_fortran_env, only: &
        sp => real32, & ! single precision, range: -3.40282347E+38 to 3.40282347E+38, smallest positive: 1.17549435E-38
        dp => real64, & ! double precision, range: -1.7976931348623157E+308 to 1.7976931348623157E+308, smallest positive: 2.2250738585072014E-308
        qp => real128, &
        i1 => int8, &   ! 1 byte integer, range: -128 to 127
        i2 => int16, &  ! 2 byte integer, range: -32768 to 32767
        i4 => int32, &  ! 4 byte integer, range: -2147483648 to 2147483647
        i8 => int64     ! 8 byte integer, range: -9223372036854775808 to 9223372036854775807
    implicit none
    private

    !> Interface for the swap subroutine.
    interface swap
        module procedure swap_int, swap_real
    end interface

    character(len=*), parameter :: fmt_general = '(*(g0,x))'
    character(len=*), parameter :: fmt_comma = '(*(g0,","))'
    character(len=*), parameter :: fmt_comma_pair = '(g0,",",g0)'

    public :: sp, dp, i1, i2, i4, i8
    public :: fmt_general, fmt_comma, count_integers_from_string
    public :: swap

contains

    subroutine swap_int(a,b)
        integer(kind=i4), intent(inout) :: a, b
        integer(kind=i4) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_int

    subroutine swap_real(a,b)
        real(kind=dp), intent(inout) :: a, b
        real(kind=dp) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_real

    function count_integers_from_string(str) result(count)
        character(len=*), intent(in) :: str
        integer(kind=i4) :: count
        integer(kind=i4) :: i, len_str
        logical :: in_number

        len_str = len_trim(str)
        count = 0
        in_number = .false.

        do i = 1, len_str
            if ((str(i:i) == ' ') .or. (str(i:i) == ',')) then
                if (in_number) then
                    count = count + 1
                    in_number = .false.
                end if
            else
                in_number = .true.
            end if
        end do

        if (in_number) count = count + 1

    end function count_integers_from_string

end module hyperSIS_kinds_mod
