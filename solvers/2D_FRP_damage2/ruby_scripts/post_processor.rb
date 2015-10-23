require 'csv'

# config

N_MODELS = 10
MODEL_PATH = '../../../models/gfrp_damage'

# module
module InpDecoder
  class << self
    def decode(inp_string)
      data = Hash.new
      lines = inp_string.each_line.map(&:strip).to_enum

      begin
        # 飛ばす
        lines.next
        lines.next

        # データ名
        data[:data_name] = lines.next

        # 節点数、要素数
        n_nodes, n_elems = lines.next.split(' ').map(&:to_i)

        # 節点
        nodes = Hash.new
        n_nodes.times do |i|
          node = lines.next.split(' ')
          nodes[node[0].to_i] = node[1..-1].map(&:to_f)
        end

        # 節点ラベル
        nodes_labels = [:x,:y,:z]

        # 要素
        elems = Hash.new
        n_elems.times do |i|
          elem = lines.next.split(' ')
          elems[elem[0].to_i] = []
          elems[elem[0].to_i] << elem[1].to_i
          elems[elem[0].to_i] << elem[2]
          elems[elem[0].to_i] += elem[3..-1].map(&:to_i)
        end

        # 要素ラベル
        elems_labels = [:material_no,:elem_type,:node1,:node2,:node3,:node4]

        # 節点データ数、要素データ数
        n_nodes_data, n_elems_data = lines.next.split(' ').map(&:to_i)
        
        # 飛ばす
        lines.next

        # 節点データラベル
        nodes_data_labels = []
        n_nodes_data.times do |i|
          nodes_data_labels << lines.next.gsub(/,*$/,"").gsub(/\s+/,"_").to_sym
        end

        # 節点データ
        headers = [:node_no] + nodes_labels + nodes_data_labels
        nodes_data = Hash.new
        n_nodes.times do |i|
          node_data = lines.next.split(' ')
          nodes_data[node_data[0].to_i] = node_data[1..-1].map(&:to_f)
        end
        data[:nodes] = CSV::Table.new(nodes.map { |i,node| CSV::Row.new(headers, [i.to_s]+node+nodes_data[i]) } )

        # 飛ばす
        lines.next

        # 要素データラベル
        elems_data_labels = []
        n_elems_data.times do |i|
          elems_data_labels << lines.next.gsub(/,*$/,"").gsub(/\s+/,"_").to_sym
        end

        # 要素データ
        headers = [:elem_no] + elems_labels + elems_data_labels
        elems_data = Hash.new
        n_elems.times do |i|
          elem_data = lines.next.split(' ')
          elems_data[elem_data[0].to_i] = elem_data[1..-1].map(&:to_f)
        end
        data[:elems] = CSV::Table.new(elems.map { |i,elem| CSV::Row.new(headers, [i.to_s]+elem+elems_data[i]) } )

        data
      rescue StopIteration
        puts "iteration reached at end"
        nil
      end
    end

  end
end

class CSV
  class Table
    def no(i)
      index = self[self.headers[0]].map(&:to_i).index(i.to_i)
      index ? self[index] : nil
    end
    def find(**conditions)
      selected = self.select do |row|
        f = true
        conditions.each do |col,val|
          f &&= (row[col] == val)
        end
        f
      end
      Table.new(selected)
    end
  end
  class Row
    def sum(h1,h2)
      i1 = self.headers.index(h1)
      i2 = self.headers.index(h2)
      self.to_a[i1..i2].inject(0) {|s,i| s += i[1].to_f }
    end
    def exp
      i1 = self.headers.index(h1)
      i2 = self.headers.index(h2)
      self.to_a[i1..i2].to_a.sum/self.length.to_f
    end
    def var
      i1 = self.headers.index(h1)
      i2 = self.headers.index(h2)
      m = self.exp
      sum = self.to_a[i1..i2].inject(0){|accum, i| accum +(i.to_f-m)**2 }
      sum/(self.length - 1).to_f
    end
    def sd
      return Math.sqrt(self.var)
    end
  end
end

module Enumerable
  def sum
    self.inject(0){|accum, i| accum + i.to_f }
  end
  def exp
    self.sum/self.length.to_f
  end
  def var
    m = self.exp
    sum = self.inject(0){|accum, i| accum +(i.to_f-m)**2 }
    sum/(self.length - 1).to_f
  end
  def sd
    return Math.sqrt(self.var)
  end
end 

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



