#
#   p_d_damaged.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   各モデルの損傷時の荷重、変位を出力するプログラム
#

require './classes.rb'
require './config.rb'

# config

#DU = 8.89e-6*2 # 変位増分
#N_MODELS = 2000 # モデル数
#MODEL_PATH = '../../../models/gfrp_damage'

# main
damage_ratios = CSV.table("#{MODEL_PATH}/damage_ratios.csv",headers: true)

io = File.open("#{MODEL_PATH}/p_d_damaged.csv","w")
io.puts "u_damaged,f_damaged"

N_MODELS.times do |i|
  reaction_forces = CSV.table("#{MODEL_PATH}/model#{i+1}/model#{i+1}.csv",headers: true)

  reaction_force_left = reaction_forces[0][0].gsub('D','e').to_f # 1ステップ目の左端荷重

  io.puts "#{damage_ratios[i][:damaged_u]}, #{reaction_force_left * damage_ratios[i][:damaged_u] / DU }"
  
end

io.close
