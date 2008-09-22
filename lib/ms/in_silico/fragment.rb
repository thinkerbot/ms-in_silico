require 'ms/in_silico/fragment_spectrum'

module Ms
  module InSilico
    # Ms::InSilico::Fragment::manifest calculates a theoretical ms/ms spectrum
    #
    # Calculates a theoretical ms/ms spectrum from a peptide sequence.
    # Configurations allow the specification of one or more fragmentation 
    # series to include, as well as charge, and intensity.  The resulting
    # data is not sorted.
    #
    #   % rap predict TVQQEL --+ dump --no-audit
    #   # date: 2008-09-15 14:37:55
    #   ---
    #   ms/in_silico/predict (26934330):
    #   - - 717.377745628191
    #     - 616.330067154091
    #     - 517.261653237891
    #     - 389.203075726491
    #     - 261.144498215091
    #     - 132.101905118891
    #     - 102.054954926291
    #     - 201.123368842491
    #     - 329.181946353891
    #     - 457.240523865291
    #     - 586.283116961491
    #     - 699.367180941891
    #
    class Fragment < Tap::Task
    
      config :series, ['y', 'b'], &c.array   # a list of the series to include
      config :charge, 1, &c.integer          # the charge for the m/z values
      config :intensity, nil, &c.num_or_nil  # a uniform intensity value
      
      def fragment_spectrum(peptide)
        FragmentSpectrum.new(peptide)
      end
      
      def process(peptide)
        spec = fragment_spectrum(peptide)
        
        masses = []
        series.each do |s| 
          masses.concat(spec.series(s))
        end
        
        masses.collect do |m|
          intensity ? [m/charge, intensity] : (m/charge)
        end
      end
    end 
  end
end