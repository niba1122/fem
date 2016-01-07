#
#   vf_at_damaged_area.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   初期損傷部位のVolume fractionをparams.csvから取得
#

require './classes.rb'
require './config.rb'


elem_no = ARGV[0].to_i || 1 # 要素を引数で指定


# 適当な1ファイル読み込み
io = File.open("#{MODEL_PATH}/model1/step1.inp","r")
data = InpDecoder.decode(io.read)
io.close

p data[:elems][7][:material_no]

SHIFT2WEFT_INDEX = []
SHIFT2WEFT_INDEX[0..13] = [0,1,2,3,4,5,6,7,8,9,10,11,12,12].map { |i| i*2 }
SHIFT2WEFT_INDEX[14..26] = [12,13,14,15,16,17,18,19,20,21,22,23,24].map { |i| i*2 }
SHIFT2WEFT_INDEX[27..51] = SHIFT2WEFT_INDEX[1..25].map { |i| i+48 }

N_PERIODS_X = 2
N_ELS_1PERIOD_X_WEFT = 96
N_ELS_1PERIOD = 704
N_ELS_1LAMINA = N_ELS_1PERIOD*N_PERIODS_X

offset_weft = LAMINATION_MISALIGNMENT.map do |shift|
  N_ELS_1PERIOD_X_WEFT/2 - ((N_ELS_1PERIOD_X_WEFT/4+SHIFT2WEFT_INDEX[shift]) % (N_ELS_1PERIOD_X_WEFT/2))
end

p LAMINATION_MISALIGNMENT
p offset_weft


LAMINATION_MISALIGNMENT.each_with_index do |shift,i|
  cnt_weft = 1
  block_no = 0
  [*1..N_ELS_1LAMINA].each do |j|
    if (data[:elems][i*N_ELS_1LAMINA+j][:material_no] == 3)
      param_no = ((cnt_weft-offset_weft[i] + N_ELS_1PERIOD_X_WEFT/2-1)/(N_ELS_1PERIOD_X_WEFT/2))+1
      if (param_no > 2*N_PERIODS_X)
        param_no = 1
      end

      mat_no_1weft = (cnt_weft + N_ELS_1PERIOD_X_WEFT/2 - offset_weft[i]) % (N_ELS_1PERIOD_X_WEFT/2)
      if (mat_no_1weft == 0)
        mat_no_1weft = N_ELS_1PERIOD_X_WEFT/2
      end

      vfmicro_no = ((mat_no_1weft-1)/4+1)*2
      if ((mat_no_1weft % 4 == 1) || (mat_no_1weft % 4 == 3))
        vfmicro_no = vfmicro_no - 1
      end

      if i*N_ELS_1LAMINA+j+1 == elem_no
        print [i+1,param_no,vfmicro_no].join(',')
      end
      
#      !material_nos((i-1)*n_els_1lamina+j) = 20+vfmicro_no
      #material_nos((i-1)*n_els_1lamina+j) = mat_no_weft(vfmicro_no,param_no,i)
#      !material_nos((i-1)*n_els_1lamina+j) = 20+param_no
      
      cnt_weft = cnt_weft + 1
    end
  end 
end
