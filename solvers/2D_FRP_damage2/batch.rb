require 'open3'

N_MODEL = 8

print "Generating random parameters"
`ruby ./random_params.rb #{N_MODEL}`

N_MODEL.times do |i|
  Open3.popen3("./2d_frp_damage2.out #{i+1}") do |stdin, stdout, stderr, w|
    stdin.close
    stdout.each do |line| print line end
    stderr.each do |line| print line end
  end
end

