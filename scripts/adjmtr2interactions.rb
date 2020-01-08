#!/usr/bin/env ruby

require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'

class OptionParse
	#
	# Return a structure describing the options.
	#
	def self.parse(args)
		# The options specified on the command line will be collected in *options*.
		# We set default values here.
		options = OpenStruct.new
		
		opts = OptionParser.new do |opts|
			opts.banner = "Usage: example.rb [options]"
		
			opts.separator ""
			opts.separator "Specific options:"
		
			# Mandatory argument.
			opts.on("-a", "--adjmtr FILE",
				"") do |adj_file_name|
				options.adj_file_name = adj_file_name
			end
			
			opts.on("-r", "--rowid FILE",
				"") do |rowid_file_name|
				options.rowid_file_name = rowid_file_name
			end

			opts.on("-c", "--columnid FILE",
				"") do |colid_file_name|
				options.colid_file_name = colid_file_name
			end

		end
		
		opts.parse!(args)
		options
	end  # parse()

end  # class OptparseExample

class KeyVal
	def initialize(key,val)
		@key=key
		@val=val
	end
	def key
		@key
	end
	def val
		@val
	end
	def report
		puts [@key,@val].join("\t")
	end
end

options = OptionParse.parse(ARGV)


rowid_file = File.open(options.rowid_file_name)
row_list = Array.new
linecount=0
rowid_file.each_with_index do |line,i|
	linecount = linecount + 1
	row_list[i] = line.strip
end

colid_file = File.open(options.colid_file_name)
col_list = Array.new
linecount=0
colid_file.each_with_index do |line,i|
	linecount = linecount + 1
	col_list[i] = line.strip
end


adj_file = File.open(options.adj_file_name)

adj_file.each_with_index do |line,i|
	line.strip().split("\t").each_with_index do |value,j|
		puts row_list[i] + "\t" + col_list[j] + "\t" + value.strip
	end
end

