#
#   damaged_position_distribution.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   指定されたモデルの初期損傷箇所分布をMicroAVSで可視化するためのファイルを生成
#

require './classes.rb'

# config

MODEL_PATH = '../../../models/gfrp_damage'
N_MODEL = 1000

# main

def get_initially_damaged_region_sigtweft(data)
  # sigT(weft)
  max_val = 0.0
  max_i = 0
  sig_t = 0.0
  min_sig_t = 1e30
  max_sig_t = -1e30
  max_sig_t_no = 0
  min_sig_t_no = 0
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

      if row[:sigII]<min_sig_t
        min_sig_t = row[:sigII]
        min_sig_t_no = row[:elem_no]
      end
      if row[:sigII]>max_sig_t
        max_sig_t = row[:sigII]
        max_sig_t_no = row[:elem_no]
      end
    end
  end
  damage_ratio_T = max_val

  #puts "## stress T (weft) ##" 
  #puts "stress at damaged element = #{sig_t} at #{max_i}"
  #puts "damage ratio = #{max_val}"
  #puts "stress_T max = #{max_sig_t} at #{max_sig_t_no}"
  #puts "stress_T min = #{min_sig_t} at #{min_sig_t_no}"
  #puts "strength_T  = #{s_t}"
  #puts "Vfmicro = #{vfmicro}"

  #puts ""

  max_i

end

def get_initially_damaged_region_sigtzweft(data)
  # sigTZ(weft)
  max_val = 0.0
  max_i = 0
  sig_tz = 0.0
  min_sig_tz = 1e30
  max_sig_tz = -1e30
  min_sig_tz_no = 0
  max_sig_tz_no = 0
  s_tz = 0.0
  vfmicro = 0.0
  damage_ratios = data[:elems].map do |row|
    if row[:material_no] == 3
      if (row[:sigII_III]/row[:max_sig_TZ]).abs > max_val
        sig_tz = row[:sigII_III]
        s_tz = row[:max_sig_TZ]
        max_val = (row[:sigII_III]/row[:max_sig_TZ]).abs
        max_i = row[:elem_no]
        vfmicro = row[:Vf]
      end

      if row[:sigII_III]<min_sig_tz
        min_sig_tz = row[:sigII_III]
        min_sig_tz_no = row[:elem_no]
      end
      if row[:sigII_III]>max_sig_tz
        max_sig_tz = row[:sigII_III]
        max_sig_tz_no = row[:elem_no]
      end
    end
  end
  damage_ratio_TZ = max_val
  #puts "## stress TZ (weft) ##"
  #puts "stress at damaged element = #{sig_tz} at #{max_i}"
  #puts "damage ratio = #{max_val}"
  #puts "stress_TZ max = #{max_sig_tz} at #{max_sig_tz_no}"
  #puts "stress_TZ min = #{min_sig_tz} at #{min_sig_tz_no}"
  #puts "strength_T  = #{s_tz}"
  #puts "Vfmicro = #{vfmicro}"

  #puts ""
  max_i

end

# あらかじめ1モデル読み込み出力用にする
io = File.open("#{MODEL_PATH}/model1/step1.inp","r")
data = InpDecoder.decode(io.read)
io.close

# 不要なデータを削除
data[:elems].headers.each do |h|
  unless [:elem_no, :material_no, :elem_type, :nodes].include? h
    data[:elems].delete h
  end
end

#data[:nodes].headers.each do |h|
#  unless [:node_no, :x, :y, :z].include? h
#    data[:nodes].delete h
#  end
#end

# メイン処理

damage_distribution = Array.new(data[:elems].length,0)
p damage_distribution.length

[*1..N_MODEL].each do |i|
  puts "Reading model#{i} ..."
  io = File.open("#{MODEL_PATH}/model#{i}/step1.inp","r")

  initially_damaged_elem_no = get_initially_damaged_region_sigtweft(InpDecoder.decode(io.read)).to_i
  damage_distribution[initially_damaged_elem_no-1] += 1

  io.close
end

data[:elems][:damage_distribution] = damage_distribution

puts "Outputting inp file..."
io = File.open("#{MODEL_PATH}/damage_distribution.inp","w")
  io.puts InpDecoder.encode(data)
io.close
