require 'csv'
require './classes.rb'

class SLS
  include Math
  attr_accessor :convergence_step
  
  N_TIMES_CONVERGENCE = 3

  def initialize(data,n_models_1step=100,max_n=10000)
    @n_models_1step = n_models_1step
    @data = data
    @max_n = max_n
    @dexp = {}
    @exp = [0.0]
    @convergence_step = nil

    threshold = @data[0...@n_models_1step].sd/(@data[0...@n_models_1step].exp*sqrt(@max_n))
    p @n_models_1step
    p @max_n
    p threshold
    p @data.length

    n_step = @data.length/@n_models_1step

    count = 0
    [*1..n_step].each do |i|
      @exp << @data[0...i*@n_models_1step].exp
      @dexp[i] = (@exp[i] - @exp[i-1]).abs

      if @dexp[i] < threshold
        count += 1
      end
      puts "#{@dexp[i]}, #{threshold}, #{count}"
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

# main
io = File.open(DATA_PATH)

data = CSV.new(io,headers: true).read

sls = SLS.new(data['damage_ratioT'].to_a,10)
p sls.convergence_step

