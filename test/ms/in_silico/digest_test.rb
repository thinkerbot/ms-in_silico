require File.join(File.dirname(__FILE__), '../../tap_test_helper.rb') 
require 'ms/in_silico/digest'

class Ms::InSilico::DigestTest < Test::Unit::TestCase
  include Ms::InSilico
  acts_as_script_test
  
  attr_accessor :d
  
  def setup
    super
    @d = Digest.new
  end
  
  def test_digest_documentation
    script_test(File.dirname(__FILE__) +  "../../../../") do |cmd|
      cmd.check "documentation", %q{
% tap run -- digest MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG --:i dump
  I[:...:]             digest MIVIGRSIVHP... to 3 peptides
MIVIGR
SIVHPYITNEYEPFAAEK
QQILSIMAG
}
    end
  end
  
  #
  # process test
  #
  
  def test_process_returns_array_of_peptide_fragments
    assert_equal %w{
      MIVIGR
      SIVHPYITNEYEPFAAEK
      QQILSIMAG
    }, d.process("MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG")
  end
  
  def test_process_skips_header_of_fasta_entries
    assert_equal %W{
      MIVIGR
      SIVHPYITN\nEYEPFAAEK
      QQILSIMAG
    }, d.process(">header\nMIVIGRSIVHPYITN\nEYEPFAAEKQQILSIMAG")
  end
end