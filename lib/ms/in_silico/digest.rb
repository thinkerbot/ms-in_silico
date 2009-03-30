require 'ms/in_silico/digester'

module Ms
  module InSilico
    # Ms::InSilico::Digest::manifest digest a protein sequence into peptides
    # Digest a protein sequence into an array of peptides.
    #
    #   % rap digest MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG --:i dump
    #     I[14:37:55]             digest MIVIGRSIVHP... to 3 peptides
    #   MIVIGR
    #   SIVHPYITNEYEPFAAEK
    #   QQILSIMAG
    #
    class Digest < Tap::Task
    
      config :digester, 'Trypsin'              # the name of the digester
      config :max_misses, 0, &c.integer        # the max # of missed cleavage sites
      config :site_digest, false, &c.boolean   # digest to sites (rather than sequences)
      
      def process(sequence)
        unless d = Digester[digester]
          raise ArgumentError, "unknown digester: #{digester}" 
        end
        
        sequence = $1 if sequence =~ /\A>.*?\n(.*)\z/m
        peptides = site_digest ? d.site_digest(sequence, max_misses): d.digest(sequence, max_misses)
        log 'digest', "#{sequence[0..10]}#{sequence.length > 10 ? '...' : ''} to #{peptides.length} peptides"
        peptides
      end
      
    end 
  end
end