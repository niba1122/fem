#
#   damaged_position.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   指定されたモデルの初期損傷箇所に関するデータを取得
#

require './classes.rb'

# config

MODEL_PATH = '../../../models/gfrp_damage'

# main

model_no = ARGV[0].to_i || 1

io = File.open("#{MODEL_PATH}/model#{model_no}/step1.inp","r")

data = InpDecoder.decode(io.read)

# sigT
max_val = 0.0
max_i = 0
sig_t = 0.0
s_t = 0.0
vfmicro = 0.0
damage_ratios = data[:elems].map do |row|
  if row[:material_no] == 3
    if (row[:sigII]/row[:max_sig_T]).abs > max_val
      sig_t = row[:sigII]
      s_t = row[:max_sig_T]
      max_val = (row[:sigII]/row[:max_sig_T]).abs
      max_i = row[:elem_no]
      vfmicro = row[:Vf]
    end
  end
end
damage_ratio_T = max_val
p data[:nodes].no(data[:elems].no(4353)[:node1])
puts "#{max_i},#{sig_t},#{s_t},#{vfmicro}"

# sigTZ
max_val = 0.0
max_i = 0
sig_tz = 0.0
s_tz = 0.0
damage_ratios = data[:elems].map do |row|
  if row[:material_no] == 3
    if (row[:sigII_III]/row[:max_sig_TZ]).abs > max_val
      sig_tz = row[:sigII_III]
      s_tz = row[:max_sig_TZ]
      max_val = (row[:sigII_III]/row[:max_sig_TZ]).abs
      max_i = row[:elem_no]
    end
  end
end
damage_ratio_TZ = max_val
puts "#{max_i},#{sig_tz},#{s_tz}"


puts InpDecoder.encode("test")
