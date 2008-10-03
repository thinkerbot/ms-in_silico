require 'constants/library'
require 'strscan'

module Ms
  module InSilico
  
    # Digester splits a protein sequence into peptides at sites specified
    # during initialization; in short Digester models a cleavage enzyme.
    # Digesters support missed cleavage sites, and can return either the
    # peptide strings or the cleavage sites.
    #
    # Digester includes {Constants::Library}[http://bioactive.rubyforge.org/constants/classes/Constants/Library.html],
    # allowing access to many common digesters using Digester[]:
    #
    #   trypsin = Digester['Trypsin']
    #   trypsin.digest('MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG')
    #   # => [
    #   # 'MIVIGR',
    #   # 'SIVHPYITNEYEPFAAEK',
    #   # 'QQILSIMAG']
    #
    #   trypsin.digest('MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG', 1)
    #   # => [
    #   # 'MIVIGR',
    #   # 'MIVIGRSIVHPYITNEYEPFAAEK',
    #   # 'SIVHPYITNEYEPFAAEK',
    #   # 'SIVHPYITNEYEPFAAEKQQILSIMAG',
    #   # 'QQILSIMAG'
    #   # ]
    #
    #   trypsin.site_digest('MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG', 1)
    #   # => [
    #   # [0,6],
    #   # [0,24],
    #   # [6,24],
    #   # [6,33],
    #   # [24,33]
    #   # ]
    #
    # ==== Enzymes
    #
    # Enzymes in the library were adapted from the default Mascot[http://www.matrixscience.com/]
    # enzyme list. Currently supported enzymes include:
    #
    # * Arg-C
    # * Asp-N
    # * Asp-N_ambic
    # * Chymotrypsin
    # * CNBr
    # * Lys-C
    # * Lys-C/P
    # * PepsinA
    # * Tryp-CNBr
    # * TrypChymo
    # * Trypsin/P
    # * V8-DE
    # * V8-E
    # * Trypsin
    # * V8-E+Trypsin
    # * V8-DE+Trypsin
    #
    # Several enzymes require two or more digesters, or functionality that
    # is not provided by Digester, and so remain unsupported:
    #
    # * CNBr+Trypsin
    # * Formic_acid
    # * LysC+AspN
    # * semiTrypsin
    #
    class Digester
      
      # The name of the digester
      attr_reader :name
      
      # A string of residues at which cleavage occurs
      attr_reader :cleave_str
      
      # A c-terminal resitriction residue which prevents 
      # cleavage at a potential cleavage site (optional).
      attr_reader :cterm_exception
      
      # True if cleavage occurs at the c-terminus of a 
      # cleavage residue, false if cleavage occurs at
      # the n-terminus.
      attr_reader :cterm_cleavage
      
      # a multiline whitespace regexp
      WHITESPACE = /\s*/m
    
      def initialize(name, cleave_str, cterm_exception=nil, cterm_cleavage=true)
        regexp = []
        0.upto(cleave_str.length - 1) {|i| regexp << cleave_str[i, 1] }
      
        @name = name
        @cleave_str = cleave_str
        @cleave_regexp = Regexp.new(regexp.join('|'))
        @cterm_exception = case 
        when cterm_exception == nil || cterm_exception.empty? then nil
        when cterm_exception.length == 1 then cterm_exception[0]
        else
          raise ArgumentError, "cterm exceptions must be a single residue: #{cterm_exception}"
        end
      
        @cterm_cleavage = cterm_cleavage
        @scanner = StringScanner.new('')
      end

      # Returns sites of digestion sites in sequence, as determined by
      # thecleave_regexp boundaries.  The digestion sites correspond
      # to the positions where a peptide begins and ends, such that
      # [n, (n+1) - n] corresponds to the [index, length] for peptide n.
      #
      #   d = Digester.new('Trypsin', 'KR', 'P')
      #   seq = "AARGGR"
      #   sites = d.cleavage_sites(seq)                 # => [0, 3, 6]
      #
      #   seq[sites[0], sites[0+1] - sites[0]]          # => "AAR"
      #   seq[sites[1], sites[1+1] - sites[1]]          # => "GGR"
      #
      # Trailing whitespace is included in the fragment.
      #
      #   seq = "AAR  \n  GGR"
      #   sites = d.cleavage_sites(seq)                 # => [0, 8, 11]
      #
      #   seq[sites[0], sites[0+1] - sites[0]]          # => "AAR  \n  "
      #   seq[sites[1], sites[1+1] - sites[1]]          # => "GGR"
      #
      # The digested section of sequence may be specified using offset 
      # and length.
      def cleavage_sites(seq, offset=0, length=seq.length-offset)
        adjustment = cterm_cleavage ? 0 : 1
        limit = offset + length
      
        positions = [offset]
        pos = scan(seq, offset, limit) do |pos|
          positions << pos - adjustment
        end

        # add the final position
        if pos < limit || positions.length == 1
          positions << limit
        end

        positions
      end

      # Returns digestion sites of sequence as [start_index, end_index] pairs,
      # allowing for missed cleavages.  Digestion sites are determined using
      # cleavage_sites; as in that method, the digested section of sequence
      # may be specified using offset and length.
      # 
      # Each [start_index, end_index] pair is yielded to the block, if given,
      # and the collected results are returned.
      def site_digest(seq, max_misses=0, offset=0, length=seq.length-offset) # :yields: start_index, end_index
        frag_sites = cleavage_sites(seq, offset, length)
      
        overlay(frag_sites.length, max_misses, 1) do |start_index, end_index|
          start_index = frag_sites[start_index]
          end_index = frag_sites[end_index]
    
          block_given? ? yield(start_index, end_index) : [start_index, end_index]
        end  
      end

      # Returns an array of peptides produced by digesting sequence, allowing for
      # missed cleavage sites. Digestion sites are determined using cleavage_sites; 
      # as in that method, the digested section of sequence may be specified using 
      # offset and length.
      def digest(seq, max_misses=0, offset=0, length=seq.length-offset)
        site_digest(seq, max_misses, offset, length).collect do |s, e|
          seq[s, e-s]
        end
      end

      protected
      
      # The cleavage regexp used to identify cleavage sites
      attr_reader :cleave_regexp # :nodoc:
      
      # The scanner used to digest strings.
      attr_reader :scanner # :nodoc:
    
      # Scans seq between offset and limit for the cleave_regexp, skipping whitespace
      # and being mindful of exception characters. The positions of the scanner at
      # each match are yielded to the block.      
      def scan(seq, offset, limit) # :nodoc:
        scanner.string = seq
        scanner.pos = offset

        while scanner.search_full(cleave_regexp, true, false)
          scanner.search_full(WHITESPACE, true, false)
          pos = scanner.pos
        
          # skip if the next character is the exception character
          next if cterm_exception != nil && seq[pos] == cterm_exception
        
          # break if you scanned past the upper limit
          break if pos > limit
        
          yield pos
        end
      
        scanner.pos
      end
    
      # Performs an overlap-collect algorithm providing the start and end 
      # indicies of spans skipping up to max_misses boundaries.
      def overlay(n, max_misses, offset) # :nodoc:
        results = []
        0.upto(n-1) do |start_index|
          0.upto(max_misses) do |n_miss|
            end_index = start_index + offset + n_miss
            break if end_index == n
    
            results << yield(start_index, end_index)
          end
        end
        results
      end
      
      #
      # Enzymes adapted from the default Mascot enzyme list.
      #
      
      class << self
        protected
        
        # Utility method to parse a mascot enzyme configuration
        # string into a Digester.
        def mascot_parse(str) # :nodoc:
          name, sense, cleave_str, cterm_exception, independent, semi_specific = str.split(/ *\t */)
          cterm_cleavage = case sense
          when 'C-Term' then true
          when 'N-Term' then false
          else raise ArgumentError, "unknown sense: #{sense}"
          end
          
          new(name, cleave_str, cterm_exception, cterm_cleavage)
        end
      end
      
      ARG_C =         mascot_parse('Arg-C 	C-Term 	R 	P 	 no 	 no')
      ASP_N =         mascot_parse('Asp-N 	N-Term 	BD 	  	no 	no')
      ASP_N_AMBIC =   mascot_parse('Asp-N_ambic 	N-Term 	DE 	  	no 	no')
      CHYMOTRYPSIN =  mascot_parse('Chymotrypsin 	C-Term 	FLWY 	P 	no 	no')
      CNBR =          mascot_parse('CNBr 	C-Term 	M 	  	no 	no')
      LYS_C =         mascot_parse('Lys-C 	C-Term 	K 	P 	no 	no')
      LYS_C_P =       mascot_parse('Lys-C/P 	C-Term 	K 	  	no 	no')
      PEPSIN_A =      mascot_parse('PepsinA 	C-Term 	FL 	  	no 	no')
      TRYP_CNBR =     mascot_parse('Tryp-CNBr 	C-Term 	KMR 	P 	no 	no')
      TRYP_CHYMO =    mascot_parse('TrypChymo 	C-Term 	FKLRWY 	P 	no 	no')
      TRYPSIN_P =     mascot_parse('Trypsin/P 	C-Term 	KR 	  	no 	no')
      V8_DE =         mascot_parse('V8-DE 	C-Term 	BDEZ 	P 	no 	no')
      V8_E =          mascot_parse('V8-E 	C-Term 	EZ 	P 	no 	no')
      TRYPSIN =       mascot_parse('Trypsin 	C-Term	KR 	P 	no 	no')
      V8_E_TRYPSIN =  mascot_parse('V8-E+Trypsin 	C-Term 	EKRZ 	P 	no 	no')
      V8_DE_TRYPSIN = mascot_parse('V8-DE+Trypsin 	C-Term 	BDEKRZ 	P 	no 	no')
      
      include Constants::Library  
      library.index_by_attribute :name
    end
  end
end