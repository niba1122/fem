#
#   classes.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   他のスクリプト用のクラス定義・拡張
#

require 'csv'

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

    #def encode(inp)
    #  str_inp = ''

    #  str_inp += "1\n"
    #  str_inp += "geom\n"
    #  str_inp += "

    #end
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

class NormRand
  include Math
  @@pi = Math.acos(-1)
  def initialize(seed,exp,sd,limit=nil)
    @exp = exp
    @sd = sd
    @limit = limit || 6.0*sd
    @prng = Random.new(seed)
  end

  def rand
    1000000.times do
      x = @prng.rand
      y = @prng.rand
      z =  sqrt(-2*log(x))*cos(2*@@pi*y) * @sd
      if z.abs <= @limit
        return z + @exp 
      end
    end
    rn
  end
end

class NormRandVfmicro
  include Math

  def initialize(seeds,exp,sd,limit=nil)
    @prngs = seeds.map { |seed| NormRand.new(seed,0.0,sd/sqrt(2.0),limit/2.0) }
    @exp = exp
    @n_vfmicro = seeds.length
  end 

  def rand
    partitions = @prngs.map do |prng|
      prng.rand
    end
    partitions << partitions[0]

    Array(0...@n_vfmicro).map { |i| @exp + partitions[i+1] - partitions[i] }
  end
end
