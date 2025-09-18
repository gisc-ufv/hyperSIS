module hyperSIS_network_mod
    use hyperSIS_kinds_mod
    use datastructs_mod, only: dynamical_list_t, fixed_list_t
    use datastructs_mod, only: log_unit, log_write, LOGGER_OK, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG
    implicit none
    private

    !> Interface for the constructor of the network type.
    interface network
        module procedure network_new
    end interface

    !> Interface for the constructor of the hyperedge type.
    !> It receives the order OR a list of nodes and returns a new hyperedge object.
    interface hyperedge
        module procedure hyperedge_new, hyperedge_new_from_nodes_list
    end interface

    !> Hyperedge object to represent a hyperedge in the network.
    !> It contains an ID, order m, and list of m + 1 nodes belonging to the hyperedge.
    type hyperedge_t
        integer(kind=i4) :: order = -1
        integer(kind=i4), allocatable :: nodes(:), dual_edges(:)
    end type

    !> Node object to represent a node in the network.
    !> It contains an ID, degree, and list of hyperedges that the node belongs to.
    type node_t
        integer(kind=i4) :: degree = 0
        integer(kind=i4), allocatable :: edges(:), dual_nodes(:)
    end type

    type network_props_t
        logical :: cleaned_null_edges = .false. ! .true. if all edges have at least a min_order
        logical :: removed_invalid_nodes_and_edges = .false. ! .true. if all nodes and edges are valid
        logical :: checked_topology_consistency = .false. ! .true. if the topology is consistent (nodes are connected to edges and vice versa)
        integer(kind=i4) :: min_order = -1 ! minimum order of the edges
    end type

    !> Network object to represent the entire hypergraph.
    !> It contains the number of nodes, number of edges, and arrays of nodes and edges objects.
    type :: network_t
        integer(kind=i4) :: num_nodes = 0
        integer(kind=i4) :: num_edges = 0
        integer(kind=i4) :: max_order = 0
        type(node_t), allocatable :: nodes(:)
        type(hyperedge_t), allocatable :: edges(:)
        type(network_props_t) :: props

        ! extra properties
        type(fixed_list_t), allocatable :: nodes_per_degree(:), edges_per_order(:)
    contains
        procedure :: init => network_init

        procedure :: build_edges_from_nodes => network_build_edges_from_nodes
        procedure :: build_nodes_from_edges => network_build_nodes_from_edges

        procedure :: print_nodes_and_edges => network_print_nodes_and_edges

        procedure :: clear_null_edges => network_clear_null_edges ! we can only remove the null edges (order -1)
        procedure :: check_topology_consistency => network_check_topology_consistency ! we can only check the topology consistency
        procedure :: remove_invalid_nodes_and_edges => network_remove_invalid_nodes_and_edges
        procedure :: clear_and_check_all => network_clear_and_check_all ! we can only check the topology consistency
        procedure :: reset_props => network_reset_props

        procedure :: destroy => network_destroy
        final :: network_finalizer
    end type

    public :: network, hyperedge
    public :: network_t, hyperedge_t, node_t

contains

    subroutine network_reset_props(net)
        class(network_t), intent(inout) :: net

        net%props%cleaned_null_edges = .false.
        net%props%removed_invalid_nodes_and_edges = .false.
        net%props%checked_topology_consistency = .false.
        net%props%min_order = -1

    end subroutine network_reset_props

    subroutine network_clear_and_check_all(net, min_order)
        class(network_t), intent(inout) :: net
        integer(kind=i4), intent(in), optional :: min_order

        net%max_order = maxval(net%edges(:)%order)

        if ((.not. net%props%cleaned_null_edges) .or. (min_order /= net%props%min_order)) call net%clear_null_edges(min_order)
        if (.not. net%props%removed_invalid_nodes_and_edges) call net%remove_invalid_nodes_and_edges()
        if (.not. net%props%checked_topology_consistency) call net%check_topology_consistency()

    end subroutine network_clear_and_check_all

    subroutine network_remove_invalid_nodes_and_edges(net)
        class(network_t), intent(inout) :: net

        ! It will remove the nodes with degree 0
        ! and edges with order -1
        ! It will be necessary to reindex everything

        ! Check nodes with degree = 0
        block
            integer(kind=i4), allocatable :: new_nodes_indexes(:)
            integer(kind=i4), allocatable :: new_edges_indexes(:)
            integer(kind=i4) :: new_num_nodes, new_num_edges, i, j, node_id, edge_id
            type(node_t), allocatable :: new_nodes(:)
            type(hyperedge_t), allocatable :: new_edges(:)

            allocate(new_nodes_indexes(net%num_nodes))
            allocate(new_nodes(net%num_nodes))

            new_num_nodes = 0
            new_nodes_indexes = -1
            do node_id = 1, net%num_nodes
                if (net%nodes(node_id)%degree /= 0) then
                    new_num_nodes = new_num_nodes + 1
                    new_nodes_indexes(node_id) = new_num_nodes

                    block
                        integer(kind=i4), allocatable :: new_edges_list(:)
                        integer(kind=i4) :: new_degree

                        allocate(new_edges_list(net%nodes(node_id)%degree))

                        new_degree = 0
                        do j = 1, net%nodes(node_id)%degree
                            edge_id = net%nodes(node_id)%edges(j)
                            if (net%edges(edge_id)%order /= -1) then ! valid edge
                                new_degree = new_degree + 1
                                new_edges_list(new_degree) = edge_id
                            end if
                        end do

                        deallocate(net%nodes(node_id)%edges)
                        net%nodes(node_id)%degree = new_degree
                        net%nodes(node_id)%edges = new_edges_list(1:new_degree)

                        new_nodes(new_num_nodes) = net%nodes(node_id)
                    end block
                end if
            end do

            deallocate(net%nodes)
            allocate(net%nodes(new_num_nodes))
            net%num_nodes = new_num_nodes
            net%nodes = new_nodes(1:new_num_nodes)
            deallocate(new_nodes)

            call log_write(LOG_DEBUG, 'New number of nodes:', net%num_nodes)

            ! Now, we need to reindex the nodes in each edge, and remove those with order -1
            allocate(new_edges_indexes(net%num_edges))
            allocate(new_edges(net%num_edges))

            new_num_edges = 0
            new_edges_indexes = -1
            do edge_id = 1, net%num_edges
                if (net%edges(edge_id)%order /= -1) then
                    new_num_edges = new_num_edges + 1
                    new_edges_indexes(edge_id) = new_num_edges

                    ! reindex the nodes in the edge
                    do i = 1, net%edges(edge_id)%order + 1
                        node_id = net%edges(edge_id)%nodes(i)
                        if (new_nodes_indexes(node_id) == -1) error stop 'should not happen! [node_id]'
                        net%edges(edge_id)%nodes(i) = new_nodes_indexes(node_id)
                    end do

                    new_edges(new_num_edges) = net%edges(edge_id)
                end if
            end do

            deallocate(net%edges)
            allocate(net%edges(new_num_edges))
            net%num_edges = new_num_edges
            net%edges = new_edges(1:new_num_edges)
            deallocate(new_edges)

            call log_write(LOG_DEBUG, 'New number of edges:', net%num_edges)

            ! Finally, update the edge_ids in the nodes
            do node_id = 1, net%num_nodes
                do i = 1, net%nodes(node_id)%degree
                    edge_id = net%nodes(node_id)%edges(i)
                    if (new_edges_indexes(edge_id) == -1) error stop 'should not happen! [edge_id]'
                    net%nodes(node_id)%edges(i) = new_edges_indexes(edge_id)
                end do
            end do
        end block

        net%props%checked_topology_consistency = .false.
        net%props%removed_invalid_nodes_and_edges = .true.
    end subroutine

    subroutine network_clear_null_edges(net, min_order_input)
        class(network_t), intent(inout) :: net
        integer(kind=i4), intent(in), optional :: min_order_input
        integer(kind=i4) :: min_order, i, j, edge_id, node_id, count
        integer(kind=i4), allocatable :: edges_with_smaller_order(:)

        if (present(min_order_input)) then
            min_order = min_order_input
        else
            min_order = 1 ! default is to have at least two nodes
        end if

        ! loop over the edges
        ! if the order is less than min_order, remove the edge
        ! and update the respective nodes
        count = 0
        do i = 1, net%num_edges
            if (net%edges(i)%order < min_order) then
                count = count + 1
            end if
        end do
        allocate(edges_with_smaller_order(count))
        count = 0
        do i = 1, net%num_edges
            if (net%edges(i)%order < min_order) then
                count = count + 1
                edges_with_smaller_order(count) = i
            end if
        end do

        call log_write(LOG_DEBUG, 'Number of edges with order smaller than', min_order, .false.)
        call log_write(LOG_DEBUG, 'is', size(edges_with_smaller_order))

        do i = 1, size(edges_with_smaller_order)
            edge_id = edges_with_smaller_order(i)

            ! we will loop over the nodes and delete the respective position there and shrink the list
            do j = 1, net%edges(edge_id)%order + 1
                node_id = net%edges(edge_id)%nodes(j)
                call clear_node_id_edges_list()
            end do

            ! remove edge
            net%edges(edge_id)%order = -1
            deallocate(net%edges(edge_id)%nodes)
        end do

        net%props%min_order = min_order
        net%props%cleaned_null_edges = .true.
        net%props%checked_topology_consistency = .false.

    contains

        subroutine clear_node_id_edges_list()
            integer(kind=i4) :: i, j, new_degree
            integer(kind=i4), allocatable :: new_edges(:)

            new_degree = net%nodes(node_id)%degree
            do i = 1, net%nodes(node_id)%degree
                if (net%nodes(node_id)%edges(i) == edge_id) then
                    net%nodes(node_id)%edges(i) = net%nodes(node_id)%edges(new_degree)
                    new_degree = new_degree - 1
                end if
            end do

            net%nodes(node_id)%degree = new_degree

            allocate(new_edges(new_degree))
            new_edges(1:new_degree) = net%nodes(node_id)%edges(1:new_degree)
            call move_alloc(new_edges, net%nodes(node_id)%edges)

        end subroutine

    end subroutine network_clear_null_edges

    subroutine network_check_topology_consistency(net)
        class(network_t), intent(inout) :: net
        integer(kind=i4) :: node_id, edge_id, i, j
        logical :: is_consistent

        is_consistent = .true.

        ! loop over nodes and check if they are in their respective edges
        do node_id = 1, net%num_nodes
            if (.not. allocated(net%nodes(node_id)%edges)) then
                if (.not. net%nodes(node_id)%degree == 0) then
                    call log_write(LOG_WARNING, 'Node', node_id, .false.)
                    call log_write(LOG_WARNING, 'has degree', net%nodes(node_id)%degree, .false.)
                    call log_write(LOG_WARNING, 'but no edges')
                    is_consistent = .false.
                    cycle
                end if
            end if
            do i = 1, net%nodes(node_id)%degree
                edge_id = net%nodes(node_id)%edges(i)
                ! check if the node is in the edge
                if (findloc(net%edges(edge_id)%nodes, node_id, dim=1) == 0) then
                    call log_write(LOG_WARNING, 'Node', node_id, .false.)
                    call log_write(LOG_WARNING, 'is not in edge', edge_id)
                    is_consistent = .false.
                end if
            end do
        end do

        ! loop over edges and check if they are in their respective nodes
        do edge_id = 1, net%num_edges
            if (.not. allocated(net%edges(edge_id)%nodes)) then
                if (.not. net%edges(edge_id)%order == -1) then
                    call log_write(LOG_WARNING, 'Edge', edge_id, .false.)
                    call log_write(LOG_WARNING, 'has order', net%edges(edge_id)%order, .false.)
                    call log_write(LOG_WARNING, 'but no nodes')
                    is_consistent = .false.
                    cycle
                end if
            end if
            do j = 1, net%edges(edge_id)%order + 1
                node_id = net%edges(edge_id)%nodes(j)
                ! check if the edge is in the node
                if (findloc(net%nodes(node_id)%edges, edge_id, dim=1) == 0) then
                    call log_write(LOG_WARNING, 'Edge', edge_id, .false.)
                    call log_write(LOG_WARNING, 'is not in node', node_id)
                    is_consistent = .false.
                end if
            end do
        end do

        net%props%checked_topology_consistency = is_consistent

        if (.not. is_consistent) then
            call log_write(LOG_ERROR, 'Topology is not consistent')
            error stop ''
        else
            call log_write(LOG_DEBUG, 'Topology is consistent')
        end if

    end subroutine network_check_topology_consistency

    subroutine network_init(net, num_nodes, num_edges)
        class(network_t) :: net
        integer(kind=i4), intent(in) :: num_nodes, num_edges

        net%num_nodes = num_nodes
        net%num_edges = num_edges
        allocate(net%nodes(num_nodes))
        allocate(net%edges(num_edges))

        net%nodes(:)%degree = 0
        net%edges(:)%order = -1
    end subroutine

    !> Constructor for the network type.
    !> It initializes the number of nodes and edges, and allocates memory for the nodes and edges arrays.
    !> @param net The network object to be initialized.
    !> @param num_nodes The number of nodes in the network.
    !> @param num_edges The number of edges in the network.
    !> @note This subroutine allocates memory for the nodes and edges arrays based on the provided sizes.
    !> @note The arrays are allocated with the sizes specified by num_nodes and num_edges.
    function network_new(num_nodes, num_edges) result(net)
        integer(kind=i4), intent(in) :: num_nodes, num_edges
        type(network_t) :: net

        call net%init(num_nodes, num_edges)

    end function network_new

    !> Constructor for the hyperedge type.
    !> It initializes the ID and order of the hyperedge, and allocates memory for the nodes array.
    !> @param order The order of the hyperedge.
    !> @note This subroutine allocates memory for the nodes array based on the provided size.
    elemental function hyperedge_new(order) result(edge)
        integer(kind=i4), intent(in) :: order
        type(hyperedge_t) :: edge

        edge%order = order
        allocate(edge%nodes(order + 1))
    end function hyperedge_new

    !> Constructor for the hyperedge type from a list of nodes.
    function hyperedge_new_from_nodes_list(nodes) result(edge)
        integer(kind=i4), intent(in) :: nodes(:)
        type(hyperedge_t) :: edge

        edge%order = size(nodes) - 1
        allocate(edge%nodes(size(nodes)))
        edge%nodes = nodes
    end function hyperedge_new_from_nodes_list

    !> Subroutine to build edges from nodes in the network.
    !> It iterates over each node and assigns the corresponding edges to the hyperedges.
    !> @param net The network object containing the nodes and edges.
    !> @note This subroutine modifies the edges of the network based on the nodes' connections.
    !> @note The edges are built by iterating over each node and assigning the corresponding edges to the hyperedges.
    !> @note The last_edge_index array is used to keep track of the last index for each hyperedge.
    subroutine network_build_edges_from_nodes(net)
        class(network_t) :: net
        integer(kind=i4) :: id_edge, i, j
        integer(kind=i4) :: last_edge_index(net%num_edges)

        last_edge_index = 0

        ! Loop over each node in the network
        do i = 1, net%num_nodes
            ! Loop over each hyperedge that the node belongs to
            do j = 1, net%nodes(i)%degree
                ! Add the node to the corresponding hyperedge
                id_edge = net%nodes(i)%edges(j)
                net%edges(id_edge)%nodes(last_edge_index(id_edge) + 1) = i
                last_edge_index(id_edge) = last_edge_index(id_edge) + 1
            end do
        end do

    end subroutine

    subroutine network_build_nodes_from_edges(net)
        class(network_t) :: net
        integer(kind=i4) :: node_id, edge_id, i
        integer(kind=i4) :: last_node_index(net%num_nodes)

        last_node_index = 0

        ! Loop over each hyperedge in the network
        do edge_id = 1, net%num_edges
            ! Loop over each node that belongs to the hyperedge
            do i = 1, net%edges(edge_id)%order + 1
                ! Add the hyperedge to the corresponding node
                node_id = net%edges(edge_id)%nodes(i)
                net%nodes(node_id)%edges(last_node_index(node_id) + 1) = edge_id
                last_node_index(node_id) = last_node_index(node_id) + 1
            end do
        end do

    end subroutine network_build_nodes_from_edges

    subroutine network_finalizer(net)
        type(network_t), intent(inout) :: net
        call net%destroy()
    end subroutine network_finalizer

    subroutine network_destroy(net)
        class(network_t), intent(inout) :: net
        integer(kind=i4) :: i

        ! Deallocate the edges and nodes arrays
        do i = 1, net%num_nodes
            if (allocated(net%nodes(i)%edges)) deallocate(net%nodes(i)%edges)
            if (allocated(net%nodes(i)%dual_nodes)) deallocate(net%nodes(i)%dual_nodes)
        end do
        do i = 1, net%num_edges
            if (allocated(net%edges(i)%nodes)) deallocate(net%edges(i)%nodes)
            if (allocated(net%edges(i)%dual_edges)) deallocate(net%edges(i)%dual_edges)
        end do

        if (allocated(net%nodes)) deallocate(net%nodes)
        if (allocated(net%edges)) deallocate(net%edges)

        net%num_nodes = 0
        net%num_edges = 0

    end subroutine network_destroy

    subroutine network_print_nodes_and_edges(net)
        class(network_t), intent(in) :: net
        integer(kind=i4) :: i

        ! Print each node
        write(*,fmt_general) 'Nodes:'
        do i = 1, net%num_nodes
            if (net%nodes(i)%degree > 0) then
                write(*,fmt_general) 'Node', i, 'degree:', net%nodes(i)%degree, 'edges:', net%nodes(i)%edges
            else
                write(*,fmt_general) 'Node', i, 'degree:', net%nodes(i)%degree, ' (isolated)'
            end if
        end do

        ! Print each hyperedge
        write(*,fmt_general) 'Hyperedges:'
        do i = 1, net%num_edges
            if (net%edges(i)%order >= 0) then
                write(*,fmt_general) 'Hyperedge', i, 'order:', net%edges(i)%order, '(', net%edges(i)%order +1, 'nodes )', &
                                                         'nodes:', net%edges(i)%nodes
            else
                write(*,fmt_general) 'Hyperedge', i, 'order:', net%edges(i)%order, ' (isolated)'
            end if
        end do

    end subroutine

end module
