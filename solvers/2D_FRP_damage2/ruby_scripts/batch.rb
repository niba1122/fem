require 'open3'

N_MODEL = 8

EXE_DIR = ".."

puts "Generating random parameters..."
`ruby ./random_params.rb #{N_MODEL}`

puts ''
puts "random parameters was Generated"
puts ''

N_MODEL.times do |i|
  Dir.chdir(EXE_DIR+"/") do
    Open3.popen3("./2d_frp_damage2.out #{i+1}") do |stdin, stdout, stderr, w|
      stdin.close
      stdout.each do |line| print line end
      stderr.each do |line| print line end
    end
  end
end

