require 'molecules/libraries/residue'
require 'constants/libraries/particle'
require 'ms/in_silico'

module Ms
  module InSilico

    # Spectrum calculates the theoretical ion series produced by a fragmentation
    # process such as collision induced disocciation (CID).  The formula used to
    # calculate the ion series were obtained from the {Matrix Science 
    # website}[http://www.matrixscience.com/].  Spectrum uses the 
    # {Constants}[http://bioactive.rubyforge.org/constants/] gem as the default
    # source of element and particle masses.
    # 
    #   spec = Ms::InSilico::Spectrum.new('TVQQEL')
    #   spec.series('b')
    #   # => [
    #   # 102.054954926291,
    #   # 201.123368842491,
    #   # 329.181946353891,
    #   # 457.240523865291,
    #   # 586.283116961491,
    #   # 699.367180941891]
    #
    #   spec.series('y')
    #   # => [
    #   # 717.377745628191,
    #   # 616.330067154091,
    #   # 517.261653237891,
    #   # 389.203075726491,
    #   # 261.144498215091,
    #   # 132.101905118891]
    #
    # ==== Formulae to Calculate Fragment Ion m/z values
    # 
    # <em>Copied directly from the Matrix Science {fragmentation help 
    # section}[http://www.matrixscience.com/help/fragmentation_help.html]</em>
    # 
    #   [N] is the molecular mass of the neutral N-terminal group, [C] is the
    #   molecular mass of the neutral C-terminal group, [M] is molecular mass
    #   of the neutral amino acid residues. To obtain m/z values, add or 
    #   subtract protons as required to obtain the required charge and divide
    #   by the number of charges. For example, to get a+, add 1 proton to the
    #   Mr value for a.  To get a--, subtract 2 protons from the Mr value for
    #   a and divide by 2. 
    #
    #    Ion Type  Neutral Mr
    #    a         [N]+[M]-CHO
    #    a*        a-NH3
    #    a°        a-H2O
    #    b         [N]+[M]-H
    #    b*        b-NH3
    #    b°        b-H2O
    #    c         [N]+[M]+NH2
    #    d         a - partial side chain
    #    v         y - complete side chain
    #    w         z - partial side chain
    #    x         [C]+[M]+CO-H
    #    y         [C]+[M]+H
    #    y*        y-NH3
    #    y°        y-H2O
    #    z         [C]+[M]-NH2
    #
    # ==== Use of alternate masses
    # By default a Spectrum will calculate the ion series' using the
    # monoisotopic masses for each element.  To calculate masses
    # differently, provide a block to new; each Element will be
    # passed to the block as needed, and the block should return
    # the element mass used in the calculation.
    #
    # Alternatively, a subclass can override the mass method; all
    # objects that need to be turned into a mass (nterm, cterm, 
    # a variety of molecules specified as strings, the elements,
    # ELECTRON, etc) are passed to mass to yield the value used
    # in any given calculation.
    #
    #--
    # ALL of the collections could be sped up using inline
    #++
    class Spectrum
      include Molecules
      include Molecules::Libraries
      include Constants::Libraries
      
      class << self
        
        def inherited(base)
          base.instance_variable_set(:@residues_to_locate, @residues_to_locate.dup)
        end
      
        # A string of residues located by scan.
        attr_accessor :residues_to_locate
      
        # Adds residues to residues_to_locate (these residues
        # will be located by scan).  Generally used when some
        # special fragmentation behavior occurs at specific
        # residues.  By default no residues are located.
        #
        #   class Subclass < Spectrum
        #     locate_residues "PS"
        #   end
        # 
        #   Subclass.new('RPPGFSPFR').residue_locations  
        #   # => {'P' => [1, 2, 6], 'S' => [5]}
        #
        # Calls to locate_residues are cumulative.
        def locate_residues(residues)
          @residues_to_locate += residues
        end
        
        # Scans the sequence to produce a ladder of masses and a
        # hash of (residue, locations) pairs which indicate the
        # indicies at which the residue occurs in sequence. The
        # ladder corresponds to the M values described above.
        #  
        # Returns [ladder, {residue => locations}].
        #
        # ==== Inputs
        # sequence:: a string
        # masses_by_byte:: an array of masses where the index of 
        #                  the mass is the byte of the
        #                  corresponding residue.
        # residues_to_locate:: a string of the residues to locate.
        #
        # Note: scan is an optimized utility function, but should
        # be replaced by an inline function to do the same.
        #
        def scan(sequence, masses_by_byte, residues_to_locate)
          locations = []
          residues_to_locate.each_byte {|byte| locations[byte] = []}

          mass = 0
          ladder = []
          sequence.each_byte do |byte|
            mass += masses_by_byte[byte]
            location = locations[byte]

            location << ladder.length if location
            ladder << mass
          end

          hash = {}
          0.upto(residues_to_locate.length-1) do |index|
            letter = residues_to_locate[index, 1]
            byte = letter[0]
            hash[letter] = locations[byte]
          end

          [ladder, hash]
        end
      end
      
      HYDROGEN = EmpiricalFormula.parse("H")
      HYDROXIDE = EmpiricalFormula.parse("OH")
      ELECTRON = Particle['Electron']
      
      self.residues_to_locate = ""
      
      # The peptide sequence.
      attr_reader :sequence
      
      # The n-terminal modification (default H)
      attr_reader :nterm
      
      # The c-terminal modification (default OH)
      attr_reader :cterm
      
      # An optional block used to calculate masses of molecules.
      attr_reader :block
      
      # A ladder of mass values corresponding to the
      # M values used in the fragmentation formulae.
      attr_reader :ladder
      
      # A hash of (residue, [locations]) pairs where
      # the locations are the indicies in sequence
      # at which residue occurs.
      attr_reader :residue_locations
    
      # Initializes a new Spectrum using the specified n- and c-terminal
      # modifications.  Masses will be calculated using the block, if
      # specified.  If no block is specified, then the monoisoptopic
      # masses will be used.
      def initialize(sequence, nterm=HYDROGEN, cterm=HYDROXIDE, &block) # :yields: element
        @sequence = sequence
        @nterm = nterm
        @cterm = cterm
        @block = block

        residue_masses = Residue.residue_index.collect do |residue| 
          next(0) if residue == nil
          mass(residue)
        end
        
        @ladder, @residue_locations = self.class.scan(
          sequence, 
          residue_masses, 
          self.class.residues_to_locate)

        @series_hash = {}
        @series_mask = {}
      end
      
      # Returns the mass of the parent ion for the sequence, given the charge.
      def parent_ion_mass(charge=1)
        (mass(nterm) + ladder.last + mass(cterm) + charge * proton_mass)/charge
      end
    
      # Returns the mass of a proton (ie Hydrogen minus an Electron)
      def proton_mass
        mass(HYDROGEN) - mass(ELECTRON)
      end
    
      # Retrieves the specfied series, assuming a charge of 1.  A different charge 
      # can be specified for the series by using '+' and '-'.  For example:
      #
      #   f = Spectrum.new 'RPPGFSPFR' 
      #   f.series('y') ==  f.y_series                      # => true
      #   f.series('b++') ==  f.b_series(2)                 # => true
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
    
      #   [N]+[M]-CHO
      def a_series(charge=1, mod=nil)
        get_series(:a, charge, mod) do
          delta = mass(mod) + mass(nterm) - mass('CHO') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      #   [N]+[M]-H
      def b_series(charge=1, mod=nil)
        get_series(:b, charge, mod) do
          delta = mass(mod) + mass(nterm) - mass('H') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      #   [N]+[M]+NH2
      def c_series(charge=1, mod=nil)
        get_series(:c, charge, mod) do
          delta = mass(mod) + mass(nterm) + mass('NH2') + charge * proton_mass
          nterm_series(delta, charge)
        end
      end

      #   [M]+H20
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

      #   [C]+[M]+CO-H
      def x_series(charge=1, mod=nil)
        get_series(:x, charge, mod) do 
          delta = mass(mod) + ladder.last + mass(cterm) + mass('CO - H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      #   [C]+[M]+H
      def y_series(charge=1, mod=nil)
        get_series(:y, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) + mass('H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      #   [C]+[M]-H
      def Y_series(charge=1, mod=nil)
        get_series(:Y, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) - mass('H') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      #   [C]+[M]-NH2
      def z_series(charge=1, mod=nil)
        get_series(:z, charge, mod) do
          delta = mass(mod) + ladder.last + mass(cterm) - mass('NH2') + charge * proton_mass
          cterm_series(delta, charge)
        end
      end

      #   [M]+H20
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
    
      # Calculates the mass of the molecule for a variety of input
      # types:
      #
      #   EmpiricalFormula   molecule.mass(&block)
      #   Particle           molecule.mass
      #   String             EmpiricalFormula.mass(molecule, &block)
      #   Numeric            molecule
      #   nil                0
      #
      def mass(molecule)
        
        # note that Particles will not actually make use of the
        # block, even though it is being passed to it.
        
        case molecule
        when EmpiricalFormula, Particle then molecule.mass(&block)
        when String then EmpiricalFormula.mass(molecule, &block)
        when nil then 0
        when Numeric then molecule
        else
          raise "cannot calculate mass of: #{molecule}"
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
      # If the series has not been initialized, the series will be 
      # initialized using the supplied block, and masked using the 
      # series_mask indicated by key (not "#{key}#{charge}").
      def get_series(key, charge=nil, mod=nil)
        series_hash["#{key}#{charge}#{mod}"] ||= mask(yield, key, mod)
      end
    
      # Mask the locations in the series by multiplying them by -1.  Mask
      # does NOT check to see if the location is negative or positive. 
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
        raise ArgumentError, "unknown series: #{s}"
      end
    end
  end
end