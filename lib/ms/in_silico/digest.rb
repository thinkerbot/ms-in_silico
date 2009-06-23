require 'tap/task'
require 'ms/in_silico/digester'

module Ms
  module InSilico
    # :startdoc::task digest a protein sequence into peptides
    # Digest a protein sequence into an array of peptides.
    #
    #   % rap digest MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG --:i dump
    #   MIVIGR
    #   SIVHPYITNEYEPFAAEK
    #   QQILSIMAG
    #
    class Digest < Tap::Task
    
      config :digester, 'Trypsin' do |digester|  # The name of the digester
        Digester[digester] or raise ArgumentError.new("unknown digester: #{digester}")
      end

      config :min_length, nil, &c.integer_or_nil # Minimum peptide length
      config :max_length, nil, &c.integer_or_nil # Maximum peptide length
      config :max_misses, 0, &c.integer          # The max # of missed cleavage sites
      config :site_digest, false, &c.boolean     # Digest to sites (rather than sequences)

      def process(sequence)
        
        # extract sequence from FASTA entries
        if sequence =~ /\A>.*?\n(.*)\z/m
          sequence = $1 
          sequence.gsub!(/\s/, "")
        end
        
        peptides = if site_digest 
          digester.site_digest(sequence, max_misses)
        else
          digester.digest(sequence, max_misses)
        end
        
        # filter
        peptides.delete_if do |peptide|
          peptide.length < min_length
        end if min_length
        
        peptides.delete_if do |peptide|
          peptide.length > max_length
        end if max_length
        
        log(:digest) { "#{sequence[0..10]}#{sequence.length > 10 ? '...' : ''} to #{peptides.length} peptides" }
        
        peptides
      end
      
    end 
  end
end