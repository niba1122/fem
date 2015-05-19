module config
  implicit none
  character(1) :: slash = "/"
  character(32) :: path_model = "../../models/" ! 最後にスラッシュを必ずつける
  double precision :: pi = dacos(-1d0)
end module
