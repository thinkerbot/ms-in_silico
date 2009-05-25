require 'tap/task'
require 'ms/in_silico/spectrum'

module Ms
  module InSilico
    
    # :startdoc::task calculates a theoretical ms/ms spectrum
    #
    # Calculates the theoretical ms/ms spectrum for a peptide sequence.
    # Configurations allow the specification of one or more fragmentation series
    # to include, as well as charge, and intensity.
    #
    # In the output the sorted fragmentation data is given first, followed by
    # a hash of header data, including the parent ion mass.
    class Fragment < Tap::Task
      
      empirical_formula_block = lambda do |value|
        case value
        when Molecules::EmpiricalFormula then value
        else Molecules::EmpiricalFormula.parse(value)
        end
      end
      
      # A block to validate a config input
      # is an EmpericalFormula.
      MOLECULE = empirical_formula_block
      
      config :series, ['y', 'b'], &c.list    # a list of the series to include
      config :charge, 1, &c.integer          # the charge for the parent ion
      config :intensity, nil, &c.numeric_or_nil  # a uniform intensity value
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