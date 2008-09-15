require 'constants/library'
require 'strscan'

module Ms
  module InSilico
    class Digester
      attr_reader :name
      attr_reader :cleave_str
      attr_reader :cleave_regexp
      attr_reader :cterm_exception
      attr_reader :cterm_cleavage
      attr_reader :scanner
    
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
          raise "cterm exceptions must be a single residue: #{cterm_exception}"
        end
      
        @cterm_cleavage = cterm_cleavage
        @scanner = StringScanner.new('')
      end

      # Returns sites of fragmentation in seq[offest, length], as determined by the 
      # cleave_regexp boundaries.  The fragmentation sites correspond to the positions 
      # where the fragments begin and end, such that sites [n, (n+1) - n]  corresponds to 
      # the [index, length] for fragment n.  Trailing whitespace is included in the fragment.
      #
      #   d = Digester.new('KR', 'P')
      #   d.fragment_sites("RGGR")                           # => [0, 1, 4]
      #   d.fragment_sites("RPPGFSPFR")                # => [0, 9]
      #   d.fragment_sites("..K..RPR..K\nPK\n\s...")   # => [0, 3, 9, 17, 20] 
      #
      def fragment_sites(seq, offset=0, length=seq.length-offset)
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

      # Returns sites of fragmentation in seq[offest, length], as determined by the 
      # cleave_regexp boundaries, allowing up to max_misses missed cleavages.  The
      # fragmentation sites correspond to the positions where the fragments begin
      # and end, such that sites [n, (n+1) - n]  corresponds to the [index, length] for 
      # fragment n. Trailing whitespace is included in the fragments.  
      #
      # site_digest yields the start and finish of each fragment to the block,  if given,  
      # and returns the collected results.
      def site_digest(seq, max_misses=0, offset=0, length=seq.length-offset) # :yields: start_index, end_index
        frag_sites = fragment_sites(seq, offset, length)
      
        overlay(frag_sites.length, max_misses, 1) do |start_index, end_index|
          start_index = frag_sites[start_index]
          end_index = frag_sites[end_index]
    
          block_given? ? yield(start_index, end_index) : [start_index, end_index]
        end  
      end

      # Returns fragments of seq[offest, length], split along the cleave_regexp
      # boundaries, allowing up to max_misses missed cleavages.  Trailing whitespace 
      # is included in the fragments.  digest yields each fragment to the block,
      # if given, and returns the collected results.
      def digest(seq, max_misses=0, offset=0, length=seq.length-offset) # :yields: fragment, start_index, end_index
        site_digest(seq, max_misses, offset, length).collect do |s, e|
          fragment = seq[s, e-s]
          block_given? ? yield(fragment, s, e) : [fragment, s, e]
        end
      end

      private
    
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

      TRYPSIN = Digester.new('Trypsin', 'KR', 'P', true)
      
      include Constants::Library  
      library.index_by_attribute :name
    end
  end
end