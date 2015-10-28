require 'csv'

#class SLS
  

DATA_PATH = '../../../models/gfrp_damage/damage_ratios.csv'

# main
io = File.open(DATA_PATH)

data = CSV.new(io).read

