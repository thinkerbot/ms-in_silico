module Ms
  
  #--
  # PRIME candidates for inline
  #++
  module InSilico
    module_function

    def fragment_scan(sequence, masses, residues="")
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
end