require 'ms/in_silico/spectrum'

module Ms
  module InSilico
    
    # Ms::InSilico::Fragment::manifest calculates a theoretical ms/ms spectrum
    #
    # Calculates the theoretical ms/ms spectrum for a peptide sequence.
    # Configurations allow the specification of one or more fragmentation series
    # to include, as well as charge, and intensity.
    #
    #   % rap fragment TVQQEL --+ dump --no-audit
    #   # date: 2008-09-15 14:37:55
    #   - - 102.054954926291
    #     - 132.101905118891
    #     - 201.123368842491
    #     - 261.144498215091
    #     - 329.181946353891
    #     - 389.203075726491
    #     - 457.240523865291
    #     - 517.261653237891
    #     - 586.283116961491
    #     - 616.330067154091
    #     - 699.367180941891
    #     - 717.377745628191
    #   - :parent_ion_mass: 717.377745628191
    #     :charge: 1
    #     :series:
    #     - y
    #     - b
    #     :nterm: H
    #     :cterm: HO
    #
    # In the output the sorted fragmentation data is given first, followed by
    # a hash of header data, including the parent ion mass.
    class Fragment < Tap::Task
      
      # A block to validate a config input
      # is an EmpericalFormula.
      MOLECULE = lambda do |value|
        case value
        when Molecules::EmpiricalFormula then value
        else Molecules::EmpiricalFormula.parse(value)
        end
      end
      
      config :series, ['y', 'b'], &c.list    # a list of the series to include
      config :charge, 1, &c.integer          # the charge for the parent ion
      config :intensity, nil, &c.num_or_nil  # a uniform intensity value
      config :nterm, 'H', &MOLECULE          # the n-terminal modification
      config :cterm, 'OH', &MOLECULE         # the c-terminal modification
      config :sort, true, &c.switch          # sorts the data by mass
      config :unmask, true, &c.switch        # remove masked (negative) masses
      
      # Constructs a hash of header data pertinent to the spectrum.
      def headers(spec)
        {
          :charge => charge,
          :parent_ion_mass => spec.parent_ion_mass(charge),
          :series => series,
          :nterm => spec.nterm.to_s,
          :cterm => spec.cterm.to_s
        }
      end
      
      def process(peptide)
        log :fragment, peptide
        spec = spectrum(peptide)
        
        masses = []
        series.each {|s| masses.concat(spec.series(s)) }
        masses.delete_if {|m| m < 0 } if unmask
        masses.sort! if sort
        masses.collect! {|m| [m, intensity] } if intensity
        
        [masses, headers(spec)]
      end
      
      protected
      
      # Returns a new Spectrum used in the calculation.
      # Primarily a hook for custom spectra in subclasses.
      def spectrum(peptide)
        Spectrum.new(peptide, nterm, cterm)
      end
      
    end 
  end
end