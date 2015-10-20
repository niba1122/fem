require 'csv'

# config

N_MODELS = 20

# module
module InpDecoder
  class << self
    def test(t)
      p 'hoge'
    end
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
        data[:n_nodes],data[:n_elems] = lines.next.split(' ').map(&:to_i)

        # 節点
        data[:nodes] = Array.new
        data[:n_nodes].times do |i|
          node = lines.next.split(' ')
          data[:nodes][node[0].to_i] = node[1..-1].map(&:to_f)
        end

        # 要素
        data[:elems] = Array.new
        data[:n_elems].times do |i|
          elem = lines.next.split(' ')
          elem_no = elem[0].to_i
          data[:elems][elem_no] = Hash.new
          
          data[:elems][elem_no][:material] = elem[1].to_i
          data[:elems][elem_no][:elem_type] = elem[2]
          data[:elems][elem_no][:nodes] = elem[3..-1].map(&:to_i)
        end

        # 節点データ数、要素データ数
        data[:n_node_data],data[:n_elem_data] = lines.next.split(' ').map(&:to_i)
        
        # 飛ばす
        lines.next

        # 節点データラベル
        node_data_labels = []
        data[:n_node_data].times do |i|
          node_data_labels << lines.next.gsub(/,*$/,"").gsub(/\s+/,"_").to_sym
        end

        # 節点データ
        data[:node_data] = Hash.new
        node_data_labels.each do |label|
          data[:node_data][label] = []
        end

        data[:n_nodes].times do |i|
          node_data = lines.next.split(' ')
          node_data[0] = node_data[0].to_i
          node_data[1..-1] = node_data[1..-1].map(&:to_f)
          node_data_labels.each_with_index do |label,i|
            data[:node_data][label][node_data[0]] = node_data[i+1]
          end
        end

        # 飛ばす
        lines.next

        # 要素データラベル
        elem_data_labels = []
        data[:n_elem_data].times do |i|
          elem_data_labels << lines.next.gsub(/,*$/,"").gsub(/\s+/,"_").to_sym
        end

        # 要素データ
        data[:elem_data] = Hash.new
        elem_data_labels.each do |label|
          data[:elem_data][label] = []
        end

        data[:n_elems].times do |i|
          elem_data = lines.next.split(' ')
          elem_data[0] = elem_data[0].to_i
          elem_data[1..-1] = elem_data[1..-1].map(&:to_f)
          elem_data_labels.each_with_index do |label,i|
            data[:elem_data][label][elem_data[0]] = elem_data[i+1]
          end
        end

        data
      rescue StopIteration
        puts "iteration reached at end"
      end
    end

  end
end

class CSV
  class Table
    def model_no(i)
      self[self[self.headers[0]].index(i.to_s)]
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

#N_MODELS.times do |i|
  
io = File.open("../../models/gfrp_damage/model1/step1.inp","r")
data = InpDecoder.decode(io.read)
io.close


io = File.open("../../models/gfrp_damage/params.csv","r")
params_csv = io.read.gsub(/\A[^\n\r]*\n/,"")
params = CSV.new(params_csv,headers: true).read

params["model No."].each do |i|
  puts params.model_no(i).sum("vfweft (1 1 1)","vfweft (1 1 24)")/24.0
end
  
io.close
