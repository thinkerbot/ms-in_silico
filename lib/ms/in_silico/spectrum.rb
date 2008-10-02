require 'molecules/libraries/residue'
require 'constants/libraries/particle'
require 'ms/in_silico'

module Ms
  module InSilico

    # Spectrum calculates the theoretical masses for ions produced by
    # fragmenting a peptide sequence in a process such as CID (collision induced 
    # disocciation).  The formula used to calculate ions for the various ion 
    # series were obtained from the MatrixScience 
    # website[http://www.matrixscience.com/help/fragmentation_help.html]
    # and have been validated against sample spectra from 
    # Mascot[http://www.matrixscience.com/] as well as 
    # ProteinProspector[http://prospector.ucsf.edu/].  One exception is listed
    # under Known Issues, below.
    # 
    # == Formulae to Calculate Fragment Ion m/z values
    # 
    # *Copied directly from {http://www.matrixscience.com/help/fragmentation_help.html}[http://www.matrixscience.com/help/fragmentation_help.html]*
    #
    # [N] is the molecular mass of the neutral N-terminal group, [C] is the molecular mass of the 
    # neutral C-terminal group, [M] is molecular mass of the neutral amino acid residues. To 
    # obtain m/z values, add or subtract protons as required to obtain the required charge and 
    # divide by the number of charges. For example, to get a+, add 1 proton to the Mr value for a. 
    # To get a--, subtract 2 protons from the Mr value for a and divide by 2. 
    #
    #  Ion Type:: 	Neutral Mr
    #  a:: 	[N]+[M]-CHO
    #  a*:: 	a-NH3
    #  a°:: 	a-H2O
    #  b:: 	[N]+[M]-H
    #  b*:: 	b-NH3
    #  b°:: 	b-H2O
    #  c:: 	[N]+[M]+NH2
    #  d:: 	a - partial side chain
    #  v:: 	y - complete side chain
    #  w:: 	z - partial side chain
    #  x:: 	[C]+[M]+CO-H
    #  y:: 	[C]+[M]+H
    #  y*:: 	y-NH3
    #  y°:: 	y-H2O
    #  z:: 	[C]+[M]-NH2
    #
    #--
    # ALL of the collections could be sped up using inline
    #++
    class Spectrum
      EmpiricalFormula = Molecules::EmpiricalFormula
      Residue = Molecules::Libraries::Residue
      Particle = Constants::Libraries::Particle
      
      class << self
        attr_reader :residues_to_locate
      
        def inherited(base)
          base.instance_variable_set(:@residues_to_locate, @residues_to_locate.dup)
        end
      
        # Specifies which residues to locate when initializing a Spectrum.
        # Most useful in subclasses.
        #
        #   class Subclass < Spectrum
        #     locate_residues "PS"
        #   end
        # 
        #   Subclass.new('RPPGFSPFR').residue_locations  
        #   # => {'P' => [1, 2, 6], 'S' => [5]}
        #
        # Calls to locate_residues are cumulative.  To reset located residues, use
        # reset_locate_residues.
        def locate_residues(residues)
          @residues_to_locate += residues
        end
      
        # Resets locate_residues such that no residues will be located.
        def reset_locate_residues
          @residues_to_locate = ""
        end
        
        def scan(sequence, masses, residues=residues_to_locate)
          locations = []
          residues.each_byte {|byte| locations[byte] = []}

          mass = 0
          ladder = []
          sequence.each_byte do |byte|
            mass += masses[byte]
            location = locations[byte]

            location << ladder.length if location
            ladder << mass
          end

          hash = {}
          0.upto(residues.length-1) do |index|
            letter = residues[index, 1]
            byte = letter[0]
            hash[letter] = locations[byte]
          end

          [ladder, hash]
        end
      end
      
      HYDROGEN = EmpiricalFormula.parse("H")
      HYDROXIDE = EmpiricalFormula.parse("OH")
      ELECTRON = Particle['Electron']
      
      reset_locate_residues
      
      # The peptide sequence 
      attr_reader :sequence
      
      # The n-terminal modification (default H)
      attr_reader :nterm
      
      # The c-terminal modification (default OH)
      attr_reader :cterm
      
      # The electron mass used in the calculation of proton_mass.
      attr_reader :electron_mass
      
      # An optional block used to calculate masses of molecules.
      attr_reader :block
      
      attr_reader :residue_masses
      attr_reader :ladder
      attr_reader :residue_locations
    
      # Initializes a new Spectrum, using the specified n- and c-terminal
      # modifications.  Masses will be calculated using the block, if
      # specified.  If no block is specified, then the monoisoptopic
      # masses will be used.
      def initialize(sequence, nterm=HYDROGEN, cterm=HYDROXIDE, residue_masses=nil, electron_mass=ELECTRON.mass, &block)
        @sequence = sequence
        @nterm = nterm
        @cterm = cterm
        @electron_mass = electron_mass
        @block = block

        @residue_masses = if residue_masses == nil
          Residue.residue_index.collect do |residue| 
            next(0) if residue == nil
            mass(residue)
          end
        else
          residue_masses
        end
        
        @ladder, @residue_locations = self.class.scan(sequence, @residue_masses)

        @series_hash = {}
        @series_mask = {}
      end
      
      # Returns the mass of the parent ion for the sequence, given the charge.
      def parent_ion_mass(charge=1)
        (mass(nterm) + ladder.last + mass(cterm) + charge * proton_mass)/charge
      end
    
      # Returns the mass of a proton (ie Hydrogen minus an Electron)
      def proton_mass
        mass(HYDROGEN) - electron_mass
      end
    
      # Retrieves the specfied series, assuming a charge of 1.  A different charge 
      # can be specified for the series by using '+' and '-'.  For example:
      #
      #   f = Spectrum.new 'RPPGFSPFR' 
      #   f.series('y') ==  f.y_series                                   # => true
      #   f.series('b++') ==  f.b_series(2)                          # => true
      #   f.series('nladder-') ==  f.nladder_series(-1)     # => true
      #
      # Series raises an error if the specified charge is zero.
      def series(s)
        s = s.to_s.strip
        case s
        when /^(immonium|nladder|cladder|[abcxyYz])(\+*)(-*)(\s[\+\-\s\w\d]+)?$/
          series = $1
          plus = $2
          minus = $3
          mod = $4.to_s.gsub(/\s/, "")
        
          charge = case
          when plus.empty? && minus.empty? then 1
          when minus.empty? then plus.length
          when plus.empty? then -minus.length
          else
            charge = plus.length - minus.length
            raise ArgumentError.new("zero charge specified in series: #{s}") if charge == 0
            charge
          end

          self.send("#{series}_series", charge, mod)
        else
          handle_unknown_series(s)
        end
      end
    
      def immonium_series(charge=1, mod=nil)
        get_series(:immonium, charge, mod) do
          delta = mass(mod) - mass('CO') 
        
          previous = 0
          series = []
          ladder.each do |current| 
            series << (current - previous + delta + charge * proton_mass)/charge
            previous = current
          end
          series
        end
      end
    
      # [N]+[M]-CHO
      def a_series(charge=1, mod=nil)
        get_series(:a, charge, mod) do
          delta = mass(mod) + mass(nterm) - mass('CHO') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      # [N]+[M]-H
      def b_series(charge=1, mod=nil)
        get_series(:b, charge, mod) do
          delta = mass(mod) + mass(nterm) - mass('H') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      # [N]+[M]+NH2
      def c_series(charge=1, mod=nil)
        get_series(:c, charge, mod) do
          delta = mass(mod) + mass(nterm) + mass('NH2') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      # [M] + H20
      #--
      # Ask Peter about these as well... Currently I'm adding water to
      # cap the ends, as if a hydrolysis reaction produced the ladder,
      # then I'm adding H for charge... is this what is intended?
      # Why not cladder[0] or cladder[-1]?
      #++
      def cladder_series(charge=1, mod=nil)
        get_series(:cladder, charge, mod) do
          delta = mass(mod) +  mass('H2O') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      #  [C]+[M]+CO-H
      def x_series(charge=1, mod=nil)
        get_series(:x, charge, mod) do 
          delta = mass(mod) + ladder.last + mass(cterm) + mass('CO - H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      # [C]+[M]+H
      def y_series(charge=1, mod=nil)
        get_series(:y, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) + mass('H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      # [C]+[M]-H
      def Y_series(charge=1, mod=nil)
        get_series(:Y, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) - mass('H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      # [C]+[M]-NH2
      def z_series(charge=1, mod=nil)
        get_series(:z, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) - mass('NH2') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      # [M] + H20
      #--
      # Ask Peter about these as well... Currently I'm adding water to
      # cap the ends, as if a hydrolysis reaction produced the ladder,
      # then I'm adding H for charge... is this what is intended?
      # Why not nladder[-1]?
      #++
      def nladder_series(charge=1, mod=nil)
        get_series(:nladder, charge, mod) do
          delta = mass(mod) + ladder.last + mass('H2O') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end
    
      protected
    
      # A hash holding all calculated series for self.  Series are keyed
      # by the type and charge of the series (ex: b1, b2, y1, y2).
      attr_accessor :series_hash
    
      # A hash holding the locations of residues that need to be masked (ie 
      # multiplied by -1) in a given series.  Mask locations should be unique
      # so that a given location will not be masked twice; the method
      # mask_locations can assist in doing so.  Series masks are keyed
      # by the series type (ex: b, y).  
      attr_accessor :series_mask
    
      # Calculates the mass of the molecule, which may be an EmpiricalFormula,
      # a string, or nil (for which the mass is 0).
      def mass(molecule)
        case molecule
        when EmpiricalFormula then molecule.mass(&block)
        when nil then 0
        when Numeric then molecule
        else EmpiricalFormula.mass(molecule, &block)
        end
      end
    
      # Generates an n-terminal series (ex: a, b, or c) by adding delta
      # to each element from ladder, and dividing by charge.  Delta,
      # therefore, should ALREADY take account of the protons added
      # by charge.
      def nterm_series(delta, charge)
        ladder.collect {|m| (m + delta)/charge }
      end
    
      # Generates a c-terminal series (ex: x, y, or z) by subtracting each
      # element from ladder from delta, and dividing by charge.  Delta,
      # therefore, should ALREADY take account of the protons added
      # by charge.
      def cterm_series(delta, charge)
        series = ladder.collect {|m| (delta - m)/charge }
        series.unshift(delta/charge)
        series.pop
        series
      end
    
      # Adds the specified locations to the series mask, ensuring that the
      # specified locations will be unique within the mask.  If overwrite
      # is true, then the input locations will overwrite any existing mask
      # locations.
      def mask_locations(series, locations, overwrite=false)
        locations = locations.collect do |location|
          location < 0 ? ladder.length + location : location
        end
      
        if overwrite
          series_mask[series] = locations.uniq
        else
          (series_mask[series] ||= []).concat(locations).uniq!
        end
      end
    
      # Retrieves the series keyed by "#{key}#{charge}" in series_hash.
      # If the series has not been initialized, the series will be initialized
      # using the supplied block, and masked using the series_mask
      # indicated by key (not "#{key}#{charge}").
      def get_series(key, charge=nil, mod=nil)
        series_hash["#{key}#{charge}#{mod}"] ||= mask(yield, key, mod)
      end
    
      # Mask the locations in the series by multiplying them by -1.  Mask does
      # NOT check to see if the location is negative or positive. 
      def mask(series, key, mod)
        locations = series_mask[key]
      
        unless mod == nil
          mod_locations = series_mask["#{key}#{mod}"] 
          if mod_locations
            locations += mod_locations 
            locations.uniq!
          end
        end
      
        locations.each {|i| series[i] *= -1} unless locations == nil
        series
      end
    
      # Hook to custom-handle an unknown series from the series method.
      # By default, handle_unknown_series raises an ArgumentError.
      def handle_unknown_series(s)
        raise ArgumentError.new("unknown series: #{s}")
      end
    end
  end
end