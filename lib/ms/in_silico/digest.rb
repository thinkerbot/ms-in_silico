require 'ms/in_silico/digester'

module Ms
  module InSilico
    # Ms::InSilico::Digest::manifest digest a protein sequence into peptides
    # Digest a protein sequence into an array of peptides.
    #
    #   % rap digest MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG --: dump --no-audit
    #     I[14:37:55]             digest MIVIGRSIVHP... to 3 peptides
    #   # date: 2008-09-15 14:37:55
    #   ---
    #   ms/in_silico/digest (23483900):
    #   - - MIVIGR
    #     - SIVHPYITNEYEPFAAEK
    #     - QQILSIMAG
    #
    class Digest < Tap::Task
    
      config :digester, 'Trypsin'              # the name of the digester
      config :max_misses, 0, &c.integer        # the max # of missed cleavage sites
      config :site_digest, false, &c.boolean   # digest to sites (rather than sequences)
      
      def process(sequence)
        unless d = Digester[digester]
          raise ArgumentError, "unknown digester: #{digester}" 
        end
        
        peptides = site_digest ? d.site_digest(sequence, max_misses): d.digest(sequence, max_misses)
        log 'digest', "#{sequence[0..10]}#{sequence.length > 10 ? '...' : ''} to #{peptides.length} peptides"
        peptides
      end
      
    end 
  end
end