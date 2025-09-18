module hyperSIS_network_io_mod
    use hyperSIS_kinds_mod
    use hyperSIS_network_mod
    use datastructs_mod, only: log_unit, log_write, LOGGER_OK, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
    implicit none
    private

    public :: network_export, network_import

contains

    subroutine network_import(net, filename, skip_check_consistency)
        class(network_t), intent(inout) :: net
        character(len=*), intent(in) :: filename
        logical, optional :: skip_check_consistency
        logical :: skip_check_consistency_choice

        if (present(skip_check_consistency)) then
            skip_check_consistency_choice = skip_check_consistency
        else
            skip_check_consistency_choice = .false.
        end if

        call network_import_edgelist(net, filename)

        if (.not. skip_check_consistency_choice) then
            call net%check_topology_consistency()
        end if

    end subroutine

    function file_exists(filename) result(res)
        character(len=*), intent(in) :: filename
        logical :: res
        integer(kind=i4) :: ios

        inquire(file=filename, exist=res, iostat=ios)
        if (ios /= 0) then
            res = .false.
        end if
    end function

    subroutine network_import_edgelist(net, filename_net)
        use, intrinsic :: iso_fortran_env, only: iostat_end
        use datastructs_mod, only: fixed_list_t, new_fixed_list_pointer

        class(network_t) :: net
        character(len=*), intent(in) :: filename_net
        integer(kind=i4) :: unit_net

        type(fixed_list_t), pointer :: edge_head => null(), edge_tail => null(), edge_current => null()

        open(newunit=unit_net, file=filename_net, status='old', action='read')

        ! read the number of nodes from the net file
        block
            integer(kind=i4) :: num_nodes, num_edges
            integer(kind=i4) :: flag
            integer(kind=i4), allocatable :: nodes(:)
            character(len=:), allocatable :: str_aux
            allocate(character(len=100000) :: str_aux)

            num_nodes = 0
            num_edges = 0
            rewind(unit_net)
            do
                read(unit_net, '(a)', iostat=flag) str_aux
                if (flag == iostat_end) exit

                allocate(nodes(count_integers_from_string(str_aux)))
                read(str_aux, *) nodes

                num_nodes = max(num_nodes, maxval(nodes))
                num_edges = num_edges + 1

                allocate(edge_current)
                call move_alloc(nodes, edge_current%list)
                edge_current%next => null()

                if (num_edges == 1) then
                    edge_head => edge_current
                    edge_tail => edge_current
                    edge_current%prev => null()
                else
                    edge_tail%next => edge_current
                    edge_current%prev => edge_tail
                    edge_tail => edge_current
                end if
            end do
            close(unit_net)
            ! allocate the nodes
            net%num_nodes = num_nodes
            net%num_edges = num_edges
            allocate(net%nodes(num_nodes))
            allocate(net%edges(num_edges))
            deallocate(str_aux)
        end block


        ! now, calculate the degree of each node and hyperedge
        block
            ! save the file to not need to read from it again later
            integer(kind=i4) :: node_id, edge_id, i
            integer(kind=i4), allocatable :: nodes(:)

            net%nodes(:)%degree = 0
            net%edges(:)%order = -1

            edge_current => edge_head

            edge_id = 1
            do while (associated(edge_current))
                call move_alloc(edge_current%list, nodes)
                net%edges(edge_id)%order = size(nodes) - 1

                do i = 1, size(nodes)
                    node_id = nodes(i)
                    net%nodes(node_id)%degree = net%nodes(node_id)%degree + 1
                end do

                call move_alloc(nodes, net%edges(edge_id)%nodes)

                edge_current => edge_current%next
                edge_id = edge_id + 1
            end do

            call destroy_list(edge_head, edge_tail)

            ! now, allocate the nodes arrays
            do node_id = 1, net%num_nodes
                allocate(net%nodes(node_id)%edges(net%nodes(node_id)%degree))
            end do

            call net%build_nodes_from_edges()

        end block

    contains

        subroutine destroy_list(head, tail)
            type(fixed_list_t), pointer :: head, tail
            type(fixed_list_t), pointer :: ptr, temp

            ptr => head
            do while (associated(ptr))
                temp => ptr%next
                if (allocated(ptr%list)) deallocate(ptr%list)
                deallocate(ptr)
                ptr => temp
            end do

            nullify(head)
            nullify(tail)
        end subroutine

    end subroutine network_import_edgelist

    subroutine network_export(net, filename)
        class(network_t), intent(in) :: net
        character(len=*), intent(in) :: filename

        call network_export_edgelist(net, filename)
    end subroutine

    !> Export the network as hyperedges list (default)
    subroutine network_export_edgelist(net, filename_net)
        use, intrinsic :: iso_fortran_env, only: iostat_end

        class(network_t) :: net
        character(len=*), intent(in) :: filename_net
        integer(kind=i4) :: unit_net
        integer(kind=i4) :: edge_id
        character(len=:), allocatable :: fmt

        fmt = fmt_general ! will never use CSV

        open(newunit=unit_net, file=filename_net, status='replace', action='write')

        do edge_id = 1, net%num_edges
            write(unit_net, fmt) net%edges(edge_id)%nodes
        end do
        close(unit_net)

    end subroutine network_export_edgelist

end module