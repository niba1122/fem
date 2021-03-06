#
#   random_params.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   Vfmicroのパラメータ値をランダムに生成する
#

require './classes.rb'

# Config

N_X_PERIODS = 2
N_Y_PERIODS = 4
EXP_MESO = 0.6
#EXP_MESO = 0.5417
SD_MESO = 0
SD_MICRO = 0.173723
#SD_MICRO = 0.1737
LIMIT_MICRO = 0.16
#LIMIT_MICRO = 0.1388 # max-ave of vfmeso

N_MICROWEFT_1WEFT = 24
N_MICROWARP_1LAMINA = N_X_PERIODS*4
N_WEFT_1LAMINA = N_X_PERIODS*2

MODEL_PATH = '../../../models/gfrp_damage'
 
ALL_CSV_PATH = "#{MODEL_PATH}/params.csv"
BASE_CSV_PATH = "#{MODEL_PATH}/params"

# Classes


# main
seed = 10
N_MODELS = ARGV[0].to_i

# warp
prng_vfwarp = Array.new
N_Y_PERIODS.times do |i|
  prng_vfwarp[i] = Array.new
  N_MICROWARP_1LAMINA.times do |j|
    prng_vfwarp[i][j] = NormRand.new(seed,EXP_MESO,EXP_MESO*SD_MICRO,LIMIT_MICRO)
    seed += 1
  end
end

# weft
prng_vfweft = Array.new
N_Y_PERIODS.times do |i|
  prng_vfweft[i] = Array.new
  N_WEFT_1LAMINA.times do |j|
    prng_vfweft[i][j] = NormRandVfmicro.new(Array(seed...(seed+N_MICROWEFT_1WEFT)),EXP_MESO,EXP_MESO*SD_MICRO,LIMIT_MICRO)
    seed += N_MICROWEFT_1WEFT
  end
end

io = File.open(ALL_CSV_PATH,"w")
io.puts N_MODELS
header = "model No."

N_Y_PERIODS.times do |i|
  N_MICROWARP_1LAMINA.times do |j|
    header << ",vfwarp (#{i+1} #{j+1})"
  end
end
p N_Y_PERIODS*N_MICROWARP_1LAMINA

N_Y_PERIODS.times do |i|
  N_WEFT_1LAMINA.times do |j|
    N_MICROWEFT_1WEFT.times do |k|
      header << ",vfweft (#{i+1} #{j+1} #{k+1})"
    end
  end
end
p N_Y_PERIODS*N_WEFT_1LAMINA*N_MICROWEFT_1WEFT
p "#{N_Y_PERIODS},#{N_WEFT_1LAMINA},#{N_MICROWEFT_1WEFT}"

io.puts header

N_MODELS.times do |k|
  Dir.mkdir(BASE_CSV_PATH) unless Dir.exist? BASE_CSV_PATH

  io_each = File.open("#{BASE_CSV_PATH}/#{k+1}.csv","w")

  vfwarp = prng_vfwarp.map { |pvs| pvs.map { |pv| pv.rand } }
  vfweft = prng_vfweft.map { |pvs| pvs.map { |pv| pv.rand } }

  str_vfwarp = vfwarp.map { |pv| pv.join(',') }.join(',')
  str_vfweft = vfweft.map { |pv| pv.map { |pv| pv.join(',') }.join(',') }.join(',')

  io.puts "#{k+1},#{str_vfwarp},#{str_vfweft}"
  io_each.puts "#{str_vfwarp},#{str_vfweft}"

end
