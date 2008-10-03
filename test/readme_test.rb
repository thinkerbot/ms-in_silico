require File.join(File.dirname(__FILE__), 'tap_test_helper.rb')
require 'ms/in_silico/digester'
require 'ms/in_silico/spectrum'

class ReadMeTest < Test::Unit::TestCase
  include Ms::InSilico
  
  def test_readme_documentation
    trypsin = Digester['Trypsin']
    peptides = trypsin.digest('MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG')
    expected = [
    'MIVIGR',
    'SIVHPYITNEYEPFAAEK',
    'QQILSIMAG']
    assert_equal expected, peptides
  
    spectrum = Spectrum.new(peptides[0])
    assert_in_delta 688.417442373391, spectrum.parent_ion_mass, 10**-12
  
    expected = [
    132.047761058391,
    245.131825038791,
    344.200238954991,
    457.284302935391,
    514.305766658991,
    670.406877687091]
    #assert_equal expected, spectrum.series('b')
  end
end