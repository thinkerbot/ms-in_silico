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
    
      config :digester, 'Trypsin'                # The name of the digester
      config :min_length, 3, &c.integer_or_nil   # Minimum peptide length
      config :max_length, nil, &c.integer_or_nil # Maximum peptide length
      config :max_misses, 0, &c.integer          # The max # of missed cleavage sites
      config :site_digest, false, &c.boolean     # Digest to sites (rather than sequences)

      def process(sequence)
        unless d = Digester[digester]
          raise ArgumentError, "unknown digester: #{digester}" 
        end
        
        # extract sequence from FASTA entries
        sequence = $1 if sequence =~ /\A>.*?\n(.*)\z/m
        
        peptides = if site_digest 
          d.site_digest(sequence, max_misses)
        else
          d.digest(sequence, max_misses)
        end
        
        # filter
        peptides.delete_if do |peptide|
          peptide.length < min_length
        end if min_length
        
        peptides.delete_if do |peptide|
          peptide.length > max_length
        end if max_length
        
        log 'digest', "#{sequence[0..10]}#{sequence.length > 10 ? '...' : ''} to #{peptides.length} peptides"
        peptides
      end
      
    end 
  end
end