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

puts "## stress T (weft) ##" 
puts "stress at damaged element = #{sig_t} at #{max_i}"
puts "damage ratio = #{max_val}"
puts "stress_T max = #{max_sig_t} at #{max_sig_t_no}"
puts "stress_T min = #{min_sig_t} at #{min_sig_t_no}"
puts "strength_T  = #{s_t}"
puts "Vfmicro = #{vfmicro}"

puts ""

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
puts "## stress TZ (weft) ##"
puts "stress at damaged element = #{sig_tz} at #{max_i}"
puts "damage ratio = #{max_val}"
puts "stress_TZ max = #{max_sig_tz} at #{max_sig_tz_no}"
puts "stress_TZ min = #{min_sig_tz} at #{min_sig_tz_no}"
puts "strength_T  = #{s_tz}"
puts "Vfmicro = #{vfmicro}"

puts ""

# sigL(warp)
max_val = 0.0
max_i = 0
sig_l = 0.0
min_sig_l = 1e30
max_sig_l = -1e30
min_sig_l_no = 0
max_sig_l_no = 0
s_l = 0.0
vfmicro = 0.0
damage_ratios = data[:elems].map do |row|
  if row[:material_no] == 2
    if (row[:sigI]/row[:max_sig_L]).abs > max_val
      sig_l = row[:sigI]
      s_l = row[:max_sig_L]
      max_val = (row[:sigI]/row[:max_sig_L]).abs
      max_i = row[:elem_no]
      vfmicro = row[:Vf]
    end

    if row[:sigI]<min_sig_l
      min_sig_l = row[:sigI]
      min_sig_l_no = row[:elem_no]
    end
    if row[:sigI]>max_sig_l
      max_sig_l = row[:sigI]
      max_sig_l_no = row[:elem_no]
    end
  end
end
damage_ratio_L = max_val

puts "## stress L (warp) ##" 
puts "stress at damaged element = #{sig_l} at #{max_i}"
puts "damage ratio = #{max_val}"
puts "stress_L max = #{max_sig_l} at #{max_sig_l_no}"
puts "stress_L min = #{min_sig_l} at #{min_sig_l_no}"
puts "strength_L  = #{s_l}"
puts "Vfmicro = #{vfmicro}"

puts ""


# sigT(warp)
max_val = 0.0
max_i = 0
sig_t = 0.0
min_sig_t = 1e30
max_sig_t = -1e30
min_sig_t_no = 0
max_sig_t_no = 0
s_t = 0.0
vfmicro = 0.0
damage_ratios = data[:elems].map do |row|
  if row[:material_no] == 2
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


puts "## stress T (warp) ##" 
puts "stress at damaged element = #{sig_t} at #{max_i}"
puts "damage ratio = #{max_val}"
puts "stress_T max = #{max_sig_t} at #{max_sig_t_no}"
puts "stress_T min = #{min_sig_t} at #{min_sig_t_no}"
puts "strength_T  = #{s_t}"
puts "Vfmicro = #{vfmicro}"

puts ""

# sigLT(warp)
max_val = 0.0
max_i = 0
sig_lt = 0.0
min_sig_lt = 1e30
max_sig_lt = -1e30
min_sig_lt_no = 0
max_sig_lt_no = 0
s_lt = 0.0
vfmicro = 0.0
damage_ratios = data[:elems].map do |row|
  if row[:material_no] == 2
    if (row[:sigI_II]/row[:max_sig_LT]).abs > max_val
      sig_lt = row[:sigI_II]
      s_lt = row[:max_sig_LT]
      max_val = (row[:sigI_II]/row[:max_sig_LT]).abs
      max_i = row[:elem_no]
      vfmicro = row[:Vf]
    end

    if row[:sigI_II]<min_sig_lt
      min_sig_lt = row[:sigI_II]
      min_sig_lt_no = row[:elem_no]
    end
    if row[:sigI_II]>max_sig_lt
      max_sig_lt = row[:sigI_II]
      max_sig_lt_no = row[:elem_no]
    end
  end
end
damage_ratio_LT = max_val
puts "## stress LT (warp) ##"
puts "stress at damaged element = #{sig_lt} at #{max_i}"
puts "damage ratio = #{max_val}"
puts "stress_LT max = #{max_sig_lt} at #{max_sig_lt_no}"
puts "stress_LT min = #{min_sig_lt} at #{min_sig_lt_no}"
puts "strength_LT  = #{s_lt}"
puts "Vfmicro = #{vfmicro}"




#puts InpDecoder.encode("test")
