module octree_mod

  implicit none

  private

  public point_type
  public octree_init
  public octree_final
  public octree_build
  public octree_update
  public octree_search

  type config_type
    integer max_num_point ! Maximum point number contained in leaf node.
    integer max_depth     ! Maximum level of branch and leaf nodes
    real(8) bbox(2, 2)
  end type config_type

  ! Points should be indexed by their id.
  type point_type
    integer id
    real(8) x(2)
  end type point_type

  ! There are two kinds of nodes:
  !   1. Branch node with children;
  !   2. Leaf node without child but containing points.
  type node_type
    integer depth
    real(8) bbox(2, 2)
    integer num_point
    integer, allocatable :: point_ids(:)
    type(node_type), pointer :: parent
    type(node_type), pointer :: children(:)
  end type node_type

  type tree_type
    type(point_type), pointer :: points(:)
    type(node_type), pointer :: root_node
  end type tree_type

  type(config_type) config
  type(tree_type) tree

contains

  subroutine octree_init(boxsize, max_num_point, max_depth, bbox)

    integer, intent(in), optional :: max_num_point
    integer, intent(in), optional :: max_depth
    real(8), intent(in), optional :: bbox(2, 2)
    real(8), intent(in) :: boxsize

    config%max_num_point = merge(max_num_point, 3, present(max_num_point))
    config%max_depth = merge(max_depth, 10, present(max_depth))
    config%bbox = merge(bbox, reshape([0.0d0, boxsize, 0.0d0, boxsize], [2, 2]), present(bbox))

    if (.not. associated(tree%root_node)) allocate(tree%root_node)
    call reset_node(tree%root_node)
    tree%root_node%depth = 1
    tree%root_node%bbox = config%bbox

  end subroutine octree_init

  subroutine octree_final() bind(c)

    call clean_node(tree%root_node)
    deallocate(tree%root_node)

  end subroutine octree_final

  subroutine build_tree(num_points, boxsize) bind(c)

    integer, intent(in) :: num_points
    integer :: i
    real(8), allocatable :: r(:), theta(:)
    real(8), intent(in) :: boxsize
    type(point_type), allocatable :: points(:)

    allocate(r(num_points))
    allocate(theta(num_points))
    allocate(points(num_points))

    print *, '   Calling octree_init'
    call octree_init(boxsize)

    print *, '   Reading coordinates of ', num_points, ' points'
    open(unit=1, file='r.txt')
    do i = 1, num_points
      read(1,*) r(i)
    end do
    close(1)
    open(unit=2, file='theta.txt')
    do i = 1, num_points
      read(2,*) theta(i)
    end do
    close(2)

    do i = 1, num_points
      points(i)%id = i
      points(i)%x(1)  = r(i)*dcos(theta(i))
      points(i)%x(2)  = r(i)*dsin(theta(i))
    end do

    print *, '   Calling octree_build'
    call octree_build(points)

  end subroutine build_tree

  recursive subroutine octree_build(points, node_)

    type(point_type), intent(in), target :: points(:)
    type(node_type), intent(inout), target, optional :: node_

    type(node_type), pointer :: node
    integer i, j
    integer num_contained_point
    type(point_type), allocatable :: contained_points(:)

    if (present(node_)) then
      node => node_
    else
      tree%points => points
      node => tree%root_node
    end if

    ! Leaf node is approached.
    if (node%depth >= config%max_depth .or. size(points) <= config%max_num_point) then
      if (size(points) > size(node%point_ids)) then
        deallocate(node%point_ids)
        allocate(node%point_ids(size(points)))
      end if
      j = 1
      do i = 1, size(points)
        if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. &
            points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2)) cycle
        node%point_ids(j) = points(i)%id
        j = j + 1
      end do
      node%num_point = j - 1
      return
    end if

    ! Copy contained points into a new array.
    num_contained_point = 0
    do i = 1, size(points)
      if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. &
          points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2)) cycle
      num_contained_point = num_contained_point + 1
    end do
    allocate(contained_points(num_contained_point))
    j = 1
    do i = 1, size(points)
      if (points(i)%x(1) < node%bbox(1, 1) .or. points(i)%x(1) > node%bbox(2, 1) .or. &
          points(i)%x(2) < node%bbox(1, 2) .or. points(i)%x(2) > node%bbox(2, 2)) cycle
      contained_points(j)%id = points(i)%id
      contained_points(j)%x = points(i)%x
      j = j + 1
    end do

    if (num_contained_point == 0) return

    ! Subdivide node and run into the child nodes.
    call subdivide_node(node)
    do i = 1, 8
      call octree_build(contained_points, node%children(i))
    end do

    ! if (node%depth == 1) then
    !   call print_tree(tree%root_node)
    ! end if

  end subroutine octree_build

  subroutine octree_update(node_)

    type(node_type), intent(inout), target, optional :: node_

    type(node_type), pointer :: node

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

  end subroutine octree_update

  subroutine tree_walk(targetx, targety, id, deltax) bind(c)

    real(8), intent(in) :: targetx, targety, targetz, deltax
    integer, intent(out) :: id
    integer num_ngb, iter
    real(8) :: x(2)
    real(8) :: searchRadius, min_dist

    searchRadius = deltax

    x(1) = targetx
    x(2) = targety

    num_ngb = 0
    iter = 0
    do while (num_ngb < 1)
      iter = iter+1
      call octree_search(x, searchRadius, num_ngb, id, searchRadius, deltax)
      searchRadius = searchRadius * 2.0d0
      ! print *, 'Found ',num_ngb,' neighbors at iteration ',iter
    end do

  end subroutine tree_walk

  recursive subroutine octree_search(x, distance, num_ngb_point, min_id, min_dist, deltax, node_)

    real(8), intent(in) :: x(2)
    real(8), intent(in) :: distance, deltax
    integer, intent(inout) :: num_ngb_point
    integer, intent(out) :: min_id
    real(8), intent(inout) :: min_dist
    type(node_type), intent(in), target, optional :: node_

    type(node_type), pointer :: node
    real(8) d2, dist2, dx(2)
    integer i

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

    if (associated(node%children)) then
      ! We are at branch node.
      do i = 1, 4
        if ((x(1) + distance) > node%children(i)%bbox(1,1) .and. &
            (x(1) - distance) < node%children(i)%bbox(2,1) .and. &
            (x(2) + distance) > node%children(i)%bbox(1,2) .and. &
            (x(2) - distance) < node%children(i)%bbox(2,2)) then
          call octree_search(x, distance, num_ngb_point, min_id, min_dist, deltax, node%children(i))
        end if
      end do
    else
      if (node%num_point == 0) return
      ! We are at leaf node.
      d2 = distance * distance
      do i = 1, node%num_point
        dx(:) = x(:) - tree%points(node%point_ids(i))%x(:)
        if ( tree%points(node%point_ids(i))%x(1) < 1.0d-10 ) then
          print *, 'INVALID TREE POINT WITH ID: ',node%point_ids(i)
          print *, '    resetting distance to mean particle spacing'
          dx(1) = deltax
        end if
        dist2 = dx(1)*dx(1) + dx(2)*dx(2) ! dot_product(dx, dx)
        if (dist2 < d2) then
          num_ngb_point = num_ngb_point + 1
          if (dist2 < min_dist*min_dist) then
            min_dist = sqrt(dist2)
            min_id = node%point_ids(i)
            ! print *, 'Found neighbor id ',min_id,' at distance ',min_dist
          end if
        end if
      end do
    end if

  end subroutine octree_search

  subroutine reset_node(node)

    type(node_type), intent(inout) :: node

    node%num_point = 0
    if (.not. allocated(node%point_ids)) allocate(node%point_ids(config%max_num_point))
    nullify(node%parent)
    if (size(node%children) == 8) deallocate(node%children)
    nullify(node%children)

  end subroutine reset_node

  subroutine subdivide_node(node)

    type(node_type), intent(inout), target :: node

    integer i, j, k, l
    real(8) bbox(2, 2)

    allocate(node%children(4))
    l = 1
    do j = 1, 2
      do i = 1, 2
        call reset_node(node%children(l))
        node%children(l)%depth = node%depth + 1
        node%children(l)%parent => node
        node%children(l)%bbox(1, 1) = node%bbox(1, 1) + (i - 1) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0
        node%children(l)%bbox(2, 1) = node%bbox(2, 1) - (2 - i) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0
        node%children(l)%bbox(1, 2) = node%bbox(1, 2) + (j - 1) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0
        node%children(l)%bbox(2, 2) = node%bbox(2, 2) - (2 - j) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0
        node%children(l)%parent => node
        l = l + 1
      end do
    end do

  end subroutine subdivide_node

  recursive subroutine clean_node(node)

    type(node_type), intent(inout) :: node

    integer i

    if (associated(node%children)) then
      do i = 1, 4
        call clean_node(node%children(i))
        deallocate(node%children(i)%point_ids)
      end do
      deallocate(node%children)
    end if

  end subroutine clean_node

  subroutine print_node(node)

    type(node_type), intent(in) :: node

    write(6, "('Bounding box: ', 6F8.2)") node%bbox
    write(6, "('Depth: ', I3)") node%depth
    write(6, "('Point number: ', I3)") node%num_point
    write(6, "('Leaf?: ', L1)") .not. associated(node%children)

  end subroutine print_node

  recursive subroutine print_tree(node)

    type(node_type), intent(in) :: node

    integer i

    if (associated(node%children)) then
      write(6, "('----------------------------------------------------------------')")
      write(6, "('Branch node: ')")
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Depth: ', I3)") node%depth
      do i = 1, 8
        call print_tree(node%children(i))
      end do
    else
      if (node%num_point == 0) return
      write(6, "('Leaf node: ')")
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Depth: ', I3)") node%depth
      write(6, "('  Points:')", advance='no')
      write(6, *) (node%point_ids(i), i = 1, node%num_point)
    end if

  end subroutine print_tree

end module octree_mod
