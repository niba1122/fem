require 'csv'
require './classes.rb'

# config

N_MODELS = 100
MODEL_PATH = '../../../models/gfrp_damage'

# main

models_data = CSV::Table.new([])
models_data_headers = [:model_no,:damage_ratioT,:damage_rationTZ,:damaged_u]

N_MODELS.times do |i|
  puts "inp file of model#{i+1} reading..."
  begin
    io = File.open("#{MODEL_PATH}/model#{i+1}/step1.inp","r")
    data = InpDecoder.decode(io.read)

    # sigT
    max_val = 0.0
    max_i = 0
    damage_ratios = data[:elems].map do |row|
      if row[:material_no] == 3
        if (row[:sigII]/row[:max_sig_T]).abs > max_val
          max_val = (row[:sigII]/row[:max_sig_T]).abs
          max_i = row[:elem_no]
        end
      end
    end
    damage_ratio_T = max_val

    # sigTZ
    max_val = 0.0
    max_i = 0
    damage_ratios = data[:elems].map do |row|
      if row[:material_no] == 3
        if (row[:sigII_III]/row[:max_sig_TZ]).abs > max_val
          max_val = (row[:sigII_III]/row[:max_sig_TZ]).abs
          max_i = row[:elem_no]
        end
      end
    end
    damage_ratio_TZ = max_val

    # damaged_u
    du = 8.89*10**(-6)
    damaged_u = du/[damage_ratio_T,damage_ratio_TZ].max

    models_data.push CSV::Row.new(models_data_headers,[i+1]<<damage_ratio_T<<damage_ratio_TZ<<damaged_u)

    io.close
#  p GC.stat
  rescue SystemCallError => e
    puts %Q(class=[#{e.class}] message=[#{e.message}])
  rescue IOError => e
    puts %Q(class=[#{e.class}] message=[#{e.message}])
  end
end
p models_data.each { |md| p md }

io = File.open("#{MODEL_PATH}/damage_ratios.csv","w")
  io.puts(models_data.to_csv)
io.close


#io = File.open("../../models/gfrp_damage/params.csv","r")
#params_csv = io.read.gsub(/\A[^\n\r]*\n/,"")
#params = CSV.new(params_csv,headers: true).read
#
#io.close



