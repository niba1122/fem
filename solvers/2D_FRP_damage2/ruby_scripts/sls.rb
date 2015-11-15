#
#   sls.rb
#   
#   (c) 2015 Nobuhito Ibaraki
#
#
#   SLS法（Stepwise Limited Sampling method）の収束判定を行うプログラム
#

require 'csv'
require './classes.rb'

class SLS
  include Math
  attr_accessor :convergence_step,:dexp,:threshold
  
  N_TIMES_CONVERGENCE = 3

  def initialize(data,n_models_1step=100,max_n=10000)
    @n_models_1step = n_models_1step
    @data = data
    @max_n = max_n
    @dexp = {}
    @exp = {}
    @convergence_step = nil

    @exp[1] = @data[0...@n_models_1step].exp
    @threshold = @data[0...@n_models_1step].sd/(@data[0...@n_models_1step].exp*sqrt(@max_n))
    p @n_models_1step
    p @max_n
    p @threshold
    p @data.length

    n_step = @data.length/@n_models_1step

    count = 0
    [*2..n_step].each do |i|
      @exp[i] = @data[0...i*@n_models_1step].exp
      @dexp[i] = ((@exp[i] - @exp[i-1])/@exp[i-1]).abs

      if @dexp[i] < @threshold
        count += 1
      end
      puts "#{@dexp[i]}, #{@threshold}, #{count}"
      if count == N_TIMES_CONVERGENCE
        @convergence_step = i
        break
      end
    end

    if count >= N_TIMES_CONVERGENCE
      return true
    else 
      return false
    end

  end

  def lsz
    
  end

end
  
#main
DATA_PATH = '../../../models/gfrp_damage/damage_ratios.csv'

io = File.open(DATA_PATH)

data = CSV.new(io,headers: true).read

# damage_ratioT
puts "----- damage_ratioT -----"
sls_damage_ratioT = SLS.new(data['damage_ratioT'].to_a,100)
p sls_damage_ratioT.convergence_step*100
p sls_damage_ratioT.dexp

sls_damage_ratioT.dexp.each do |step,dexp|
  puts "#{step},#{dexp}"
end

# damage_ratioTZ
puts "----- damage_ratioTZ -----"
sls_damage_ratioTZ = SLS.new(data['damage_rationTZ'].to_a,100)
p sls_damage_ratioTZ.convergence_step*100
p sls_damage_ratioTZ.dexp

sls_damage_ratioTZ.dexp.each do |step,dexp|
  puts "#{step},#{dexp}"
end

# u at damaged
puts "----- u at damaged -----"
sls_damaged_u = SLS.new(data['damaged_u'].to_a,100)
p sls_damaged_u.convergence_step*100
p sls_damaged_u.dexp

sls_damaged_u.dexp.each do |step,dexp|
  puts "#{step},#{dexp}"
end
